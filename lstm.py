import tensorflow as tf
import numpy as np
import pprint
pp = pprint.PrettyPrinter(indent=4)
sess = tf.InteractiveSession()


# train Parameters
seq_length = 9
data_dim = 1
hidden_dim = 30
output_dim = 2
num_classes = 500
learning_rate = 0.001
iterations = 2000

xy = np.loadtxt('data_for_train2.csv', delimiter=',',dtype=np.float32)
XY = np.loadtxt('data_for_test2.csv',delimiter=',',dtype=np.float32)
np.random.shuffle(xy)
np.random.shuffle(XY)
x = xy[:,1:-1]
y = xy[:, -1:]
x_test = XY[:, 1:-1]
y_test = XY[:,-1:]

# build a dataset
dataX = []
dataY = []
ValX = []
ValY = []
testX = []
testY = []
for i in range(0, len(y)):
    _x = [[z] for z in x[i][:]]
    _y = y[i]  # snp
    print(_x, "->", _y)
    if i < (0.8*len(y)):
        dataX.append(_x)
        dataY.append(_y)
    # print(dataX)
    else:
        ValX.append(_x)
        ValY.append(_y)

for i in range(0,len(y_test)):
    _x = [[z] for z in x_test[i][:]]
    _y =  y_test[i]
    testX.append(_x)
    testY.append(_y)

#
# input place holders
X = tf.placeholder(tf.float32, [None, seq_length, data_dim])
Y = tf.placeholder(tf.int32, [None, 1])
Y_one_hot = tf.one_hot(Y, 2)
Y_one_hot = tf.reshape(Y_one_hot, [-1,2])

# build a LSTM network
cell = tf.contrib.rnn.BasicLSTMCell(
    num_units=hidden_dim, state_is_tuple=True, activation=tf.tanh)
outputs, _states = tf.nn.dynamic_rnn(cell, X, dtype=tf.float32)

FullC_output = tf.contrib.layers.fully_connected(
    outputs[:, -1], num_classes, activation_fn=None)  # We use the last cell's output
print("FullC_output",FullC_output)
weights = tf.Variable(tf.random_normal([num_classes,output_dim]))
bias = tf.Variable(tf.random_normal([output_dim]))
prediction=tf.nn.softmax(tf.matmul(FullC_output,weights)+bias)
print("prediction",prediction)
# cost/loss
# loss = tf.reduce_sum(tf.square(Y_pred - Y))  # sum of the squares
loss = tf.reduce_mean(tf.nn.softmax_cross_entropy_with_logits(logits=prediction, labels=Y_one_hot))
# optimizer
optimizer = tf.train.AdamOptimizer(learning_rate)
train = optimizer.minimize(loss)

#predict
prediction = tf.argmax(prediction,1)
correct_prediction = tf.equal(prediction, tf.argmax(Y_one_hot,1))
accuracy = tf.reduce_mean(tf.cast(correct_prediction, tf.float32))

with tf.Session() as sess:
    init = tf.global_variables_initializer()
    sess.run(init)

    # Training step
    for i in range(iterations):
        _, step_loss = sess.run([train, loss], feed_dict={
                                X: dataX, Y: dataY})
        # print(sess.run([train, loss, X.eval(), Y.eval()],feed_dict={X:dataX, Y:dataY}))
        print("[step: {}] loss: {}".format(i, step_loss))
        # val_pred = sess.run()

    # Test step
    test_predict = sess.run(prediction, feed_dict={X: testX})
    for p, y in zip(test_predict, testY):
        print("[{}] Prediction : {} True Y: {}".format(p == int(y), p, int(y)))

    print(sess.run(accuracy, feed_dict={X: testX, Y : testY}))