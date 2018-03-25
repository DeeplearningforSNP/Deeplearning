import tensorflow as tf
import numpy as np
sess=tf.Session()
hidden_size=200
batch_size=20
num_classes=200
seq_leng=9
out_dim=2

#retfile=open('val.csv','w')

xy_train=np.loadtxt('data_for_train.csv',delimiter=', ',dtype=np.float32)
np.random.shuffle(xy_train)
'''
x_train=xy_train[0:640,1:-2].reshape(-1,seq_leng,1)
y_train=xy_train[0:640,-2:]
'''
xy_val=np.loadtxt('data_for_validation.csv',delimiter=', ',dtype=np.float32)
np.random.shuffle(xy_val)
''''''
x_val=xy_val[0:300,1:-2].reshape(-1,seq_leng,1)
y_val=xy_val[0:300,-2:]
''''''
xy_test=np.loadtxt('data_for_test.csv',delimiter=', ',dtype=np.float32)
np.random.shuffle(xy_test)

x_test=xy_test[0:652,1:-2].reshape(-1,seq_leng,1)
y_test=xy_test[0:652,-2:]

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

correct_prediction=tf.equal(tf.argmax(prediction,1),tf.argmax(Y,1))
accuracy=tf.reduce_mean(tf.cast(correct_prediction,tf.float32))

ret_matrix=tf.Variable(tf.zeros(2,2))

last_loss=99999999.
sess.run(tf.global_variables_initializer())
for epoch in range(200):
    for i in range(int(900/batch_size)):
        xy_batch=xy_train[i*batch_size:i*batch_size+batch_size]
        np.random.shuffle(xy_batch)
        x_batch,y_batch = xy_batch[:,1:-2].reshape(-1,seq_leng,1),xy_batch[:,-2:]

        l,_=sess.run([loss,train],feed_dict={X:x_batch,Y:y_batch})
    if epoch%10==0:
#        retfile.write("%s %s\n"%(sess.run(accuracy,feed_dict={X:x_val,Y:y_val}),sess.run(loss,feed_dict={X:x_val,Y:y_val})))
        val_loss=sess.run(loss,feed_dict={X:x_val,Y:y_val})
        if last_loss < val_loss:
            break
        last_loss=val_loss
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
'''
#for i in range(len(y_test)):
#    retfile.write("%s\t%s\n"%(sess.run(tf.argmax(Y[i]),feed_dict={Y:y_test}),sess.run(tf.argmax(prediction[i]),feed_dict={X:x_test})))
#print(sess.run([tf.argmax(prediction,1),tf.cast(Y,tf.float32)],feed_dict={X:x_test,Y:y_test}))
print('prediction',sess.run(correct_prediction,feed_dict={X:x_test,Y:y_test}))
print('accuracy',sess.run(accuracy,feed_dict={X:x_test,Y:y_test}))
for i in range(len(y_test)):
#    print(sess.run([Y[i][0],tf.argmax(prediction,1)[i]],feed_dict={Y:y_test,X:x_test}))
    if sess.run(Y[i][0],feed_dict={Y:y_test})==1:
        if sess.run(tf.argmax(prediction,1)[i],feed_dict={X:x_test})==1:
            ret_matrix+=[[0,0],[1,0]]
        else:
            ret_matrix+=[[0,0],[0,1]]
    else:
        if sess.run(tf.argmax(prediction,1)[i],feed_dict={X:x_test})==0:
            ret_matrix+=[[0,1],[0,0]]
        else:
            ret_matrix+=[[1,0],[0,0]]
print(sess.run(ret_matrix,feed_dict={X:x_test,Y:y_test}))

#retfile.close()
