import tensorflow as tf
import numpy as np
sess=tf.Session()
hidden_size=200
batch_size=10
num_classes=200
seq_leng=9
out_dim=2

retfile=open('val.csv','w')

xy_train=np.loadtxt('data_for_train.csv',delimiter=', ',dtype=np.float32)
np.random.shuffle(xy_train)
'''
x_train=xy_train[0:900,1:-2].reshape(-1,seq_leng,1)
y_train=xy_train[0:900,-2:]
'''
xy_val=np.loadtxt('data_for_validation.csv',delimiter=', ',dtype=np.float32)
np.random.shuffle(xy_val)
''''''
x_val=xy_val[0:300,1:-2].reshape(-1,seq_leng,1)
y_val=xy_val[0:300,-2:]
''''''
xy_test=np.loadtxt('data_for_test.csv',delimiter=', ',dtype=np.float32)
np.random.shuffle(xy_test)

x_test=xy_test[0:600,1:-2].reshape(-1,seq_leng,1)
y_test=xy_test[0:600,-2:]

X=tf.placeholder(tf.float32,shape=[None,seq_leng,1])
Y=tf.placeholder(tf.int64,shape=[None,2])
learning_rate=0.001

cell = tf.contrib.rnn.BasicLSTMCell(num_units=hidden_size, state_is_tuple=True)
outputs, _states = tf.nn.dynamic_rnn(cell, X, dtype=tf.float32)
fc_out=tf.contrib.layers.fully_connected(outputs[:,-1],num_classes,activation_fn=None)

weights = tf.Variable(tf.random_normal([num_classes,out_dim]))
bias = tf.Variable(tf.random_normal([out_dim]))
prediction=tf.nn.softmax(tf.matmul(fc_out,weights)+bias)
loss=tf.reduce_mean(tf.square(tf.cast(Y,tf.float32)-prediction))
train = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(loss)

correct_prediction=tf.equal(prediction,tf.cast(Y,tf.float32))
accuracy=tf.reduce_mean(tf.cast(correct_prediction,tf.float32))

w=[]
b=[]

sess.run(tf.global_variables_initializer())
for epoch in range(2000):
    for i in range(int(900/batch_size)):
        xy_batch=xy_train[i*batch_size:i*batch_size+batch_size]
        np.random.shuffle(xy_batch)
        x_batch,y_batch = xy_batch[:,1:-2].reshape(-1,seq_leng,1),xy_batch[:,-2:]

        l,_=sess.run([loss,train],feed_dict={X:x_batch,Y:y_batch})
#        if i%3==0:
#    print('val accuracy : ',
    retfile.write("%s\n"%sess.run(accuracy,feed_dict={X:x_val,Y:y_val}))

    w.insert(epoch,weights.eval(session=sess))
    b.insert(epoch,bias.eval(session=sess))
    print('weight\n',w[epoch])
    print('bias\n',b[epoch])
    print('result : ',epoch,l)
'''
last_loss=99999999
for epoch in range(2000):
    for i in range(int(300/batch_size)):
        xy_batch=xy_val[i*batch_size:i*batch_size+batch_size]
        np.random.shuffle(xy_batch)
        x_batch,y_batch = xy_batch[:,1:-2].reshape(-1,seq_leng,1),xy_batch[:,-2:]

        l=sess.run(loss,feed_dict={X:x_batch,Y:y_batch})
    print(epoch,l)
    if last_loss*10 < l:
        break
    last_loss=l
snpx=0
xsnp=0
for i in range(len(Y_test)):
    if Y[i][0]==1 and tf.argmax(prediction,1)[i]==1:
        snpx+=1
    elif Y[i][0]==0 and tf.argmax(prediction,1)[i]==0:
        xsnp+=1
#correct_prediction=tf.equal(prediction,tf.cast(Y,tf.float32))
#accuracy=tf.reduce_mean(tf.cast(correct_prediction,tf.float32))
print(sess.run(accuracy,feed_dict={X:x_test,Y:y_test}))
sess.run(tf.global_variables_initializer())
print(sess.run([xsnp,snpx],feed_dict={X:x_test,Y:y_test}))
'''
retfile.close()
