import tensorflow as tf
import numpy as np
from collections import Counter

sess=tf.Session()
hidden_size=200
batch_size=20
num_classes=200
seq_leng=9
out_dim=2

xy_train=np.loadtxt('train_qsum.csv',delimiter=', ',dtype=np.float32)

xy_val=np.loadtxt('val_qsum.csv',delimiter=', ',dtype=np.float32)
np.random.shuffle(xy_val)
x_val=xy_val[0:2000,1:-2].reshape(-1,seq_leng,1)
y_val=xy_val[0:2000,-2:]

xy_test=np.loadtxt('test_qsum.csv',delimiter=', ',dtype=np.float32)
np.random.shuffle(xy_test)
x_test=xy_test[0:5063,1:-2].reshape(-1,seq_leng,1)
y_test=xy_test[0:5063,-2:]

X=tf.placeholder(tf.float32,shape=[None,seq_leng,1])
Y=tf.placeholder(tf.int64,shape=[None,2])
learning_rate=0.0001

cell = tf.contrib.rnn.BasicLSTMCell(num_units=hidden_size, state_is_tuple=True)
outputs, _states = tf.nn.dynamic_rnn(cell, X, dtype=tf.float32)
fc_out=tf.contrib.layers.fully_connected(outputs[:,-1],num_classes,activation_fn=tf.nn.relu)

weights = tf.get_variable("weights",shape=[num_classes,out_dim],initializer=tf.contrib.layers.xavier_initializer())
bias = tf.get_variable("bias",shape=[out_dim],initializer=tf.contrib.layers.xavier_initializer())
logits=tf.matmul(fc_out,weights)+bias
prediction=tf.nn.softmax(logits)
cost=tf.nn.softmax_cross_entropy_with_logits(logits=logits,labels=Y)
loss=tf.reduce_mean(cost)
loss_summ=tf.summary.scalar("loss",loss)
train = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(loss)

correct_prediction=tf.equal(tf.argmax(prediction,1),tf.argmax(Y,1))
accuracy=tf.reduce_mean(tf.cast(correct_prediction,tf.float32))

ret_matrix=tf.Variable(tf.zeros(2,2))

merged_summ=tf.summary.merge_all()
writer=tf.summary.FileWriter("./logs/snp_loss")
writer.add_graph(sess.graph)
sess.run(tf.global_variables_initializer())
for epoch in range(500):
    np.random.shuffle(xy_train)
    for i in range(int(18000/batch_size)):
        xy_batch=xy_train[i*batch_size:i*batch_size+batch_size]
        x_batch,y_batch = xy_batch[:,1:-2].reshape(-1,seq_leng,1),xy_batch[:,-2:]

        l,_=sess.run([loss,train],feed_dict={X:x_batch,Y:y_batch})
    print(epoch, l)
    summary=sess.run(merged_summ,feed_dict={X:x_val,Y:y_val})
    writer.add_summary(summary,global_step=epoch)

print('prediction',sess.run(correct_prediction,feed_dict={X:x_test,Y:y_test}))
print('accuracy',sess.run(accuracy,feed_dict={X:x_test,Y:y_test}))

Y_real,Y_pred = np.array(sess.run([tf.argmax(Y,1),tf.argmax(prediction,1)],feed_dict={X:x_test,Y:y_test}))
Y_ret=2*Y_real+Y_pred
print(Counter(Y_ret))
