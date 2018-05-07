import tensorflow as tf
import numpy as np
from collections import Counter

sess=tf.Session()
hidden_size=200
batch_size=100
num_classes=200
seq_leng=10
out_dim=2

xy_train=np.loadtxt('train_qsum.csv',delimiter=', ',dtype=np.float32)

xy_val=np.loadtxt('val_qsum.csv',delimiter=', ',dtype=np.float32)
np.random.shuffle(xy_val)
x_val=xy_val[0:2000,1:-2].reshape(-1,seq_leng,1)
y_val=xy_val[0:2000,-2:]

test=open('test_qsum.csv','r')
xy_test=test.read().split('\n')
pos_dic={}
test_pos=[]
test_datas=[]
for i in range(len(xy_test)-1):
    pos_dic[int(xy_test[i].split(', ')[0])]=np.array(xy_test[i].split(', ')[1:],dtype=np.float32)
    xy_test[i]=xy_test[i].split(', ')[1:]
xy_test=np.array(xy_test[:-1],dtype=np.float32)
np.random.shuffle(xy_test)

X=tf.placeholder(tf.float32,shape=[None,seq_leng,1],name='X')
Y=tf.placeholder(tf.int64,shape=[None,2])
learning_rate=0.00002

cell = tf.contrib.rnn.BasicLSTMCell(num_units=hidden_size, state_is_tuple=True)
outputs, _states = tf.nn.dynamic_rnn(cell, X, dtype=tf.float32)
fc_out=tf.contrib.layers.fully_connected(outputs[:,-1],num_classes,activation_fn=tf.nn.relu)

weights = tf.get_variable("weights",shape=[num_classes,out_dim],initializer=tf.contrib.layers.xavier_initializer())
bias = tf.get_variable("bias",shape=[out_dim],initializer=tf.contrib.layers.xavier_initializer())
logits=tf.matmul(fc_out,weights)+bias
prediction=tf.nn.softmax(logits,name='prediction')
cost=tf.nn.softmax_cross_entropy_with_logits(logits=logits,labels=Y)
loss=tf.reduce_mean(cost)
loss_summ=tf.summary.scalar("loss",loss)
train = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(loss)

correct_prediction=tf.equal(tf.argmax(prediction,1),tf.argmax(Y,1))
accuracy=tf.reduce_mean(tf.cast(correct_prediction,tf.float32))

merged_summ=tf.summary.merge_all()
writer=tf.summary.FileWriter("./logs/sum")
writer.add_graph(sess.graph)

saver=tf.train.Saver()
sess.run(tf.global_variables_initializer())

overfitting=0
min_loss=9999999
for epoch in range(100000):
    np.random.shuffle(xy_train)
    for i in range(int(18000/batch_size)):
        xy_batch=xy_train[i*batch_size:i*batch_size+batch_size]
        x_batch,y_batch = xy_batch[:,1:-2].reshape(-1,seq_leng,1),xy_batch[:,-2:]

        l,_=sess.run([loss,train],feed_dict={X:x_batch,Y:y_batch})
    print(epoch, l)
    if epoch%20==0:
        saver.save(sess,'./learning_model/sum.ckpt')

    val_loss=sess.run(loss,feed_dict={X:x_val,Y:y_val})
    if val_loss < min_loss:
        overfitting=0
        min_loss=val_loss
    else:
        overfitting+=1
    if overfitting>10:
        break

    summary=sess.run(merged_summ,feed_dict={X:x_val,Y:y_val})
    writer.add_summary(summary,global_step=epoch)

ret=open('training_ret.txt','w')
Y_real=[]
Y_pred=[]
accu=0
for i in range(int(np.ceil(5063/batch_size))):
    test_batch=xy_test[i*batch_size:i*batch_size+batch_size]
    test_x,test_y=test_batch[:,:-2].reshape(-1,seq_leng,1),test_batch[:,-2:]

    accu+=sess.run(accuracy,feed_dict={X:test_x,Y:test_y})

    Y_real+=list(sess.run(tf.argmax(Y,1),feed_dict={Y:test_y}))
    Y_pred+=list(sess.run(tf.argmax(prediction,1),feed_dict={X:test_x}))

print('accuracy : ',accu/int(np.ceil(5063/batch_size)))

Y_ret=2*np.array(Y_real)+np.array(Y_pred)
print(Counter(Y_ret))
true_snp=[i for i, x in enumerate(Y_ret) if x==0]
false_snp=[i for i, x in enumerate(Y_ret) if x==1]

for i in range(len(true_snp)):
    for pos in list(pos_dic.keys()):
        if np.array_equal(xy_test[true_snp[i]],pos_dic[pos]):
            ret.write("%s\n"%pos)
ret.write("\n\n")
for i in range(len(false_snp)):
    for pos in list(pos_dic.keys()):
        if np.array_equal(xy_test[false_snp[i]],pos_dic[pos]):
            ret.write("%s\n"%pos)
test.close()
ret.close()
