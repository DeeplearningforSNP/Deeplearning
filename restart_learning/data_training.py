import tensorflow as tf
import numpy as np
from tensorflow.contrib import layers
from collections import Counter

sess=tf.Session()
hidden_size=200
batch_size=100
num_classes=200
seq_leng=9
out_dim=2

xy_train=np.loadtxt('learn_train.csv',delimiter=', ',dtype=np.float32)

xy_val=np.loadtxt('learn_val.csv',delimiter=', ',dtype=np.float32)
np.random.shuffle(xy_val)
x_val=xy_val[0:2000,1:-2].reshape(-1,seq_leng,1)
y_val=xy_val[0:2000,-2:]

xy_test=np.loadtxt('learn_test.csv',delimiter=', ',dtype=np.float32)
np.random.shuffle(xy_test)

X=tf.placeholder(tf.float32,shape=[None,seq_leng,1],name='X')
Y=tf.placeholder(tf.int64,shape=[None,2])
learning_rate=0.00002

#   lstm으로 hidden size=256 을 many to one 방식으로 뽑아
#   fully connected에 전달하여 num classes=256 으로 진행시킴(기존 코드는 200 200)
#   fully connected에서 activation fn은 relu이고 weight와 bias 초기화는 기존 코드에서 추가됨
#   이후 weight와 bias값을 통해 fully connected를 진행하는 방식에서
#   activation fn이 없어 의미가 없다고 보고 이 부분은 fully connected로 변경(out dim=2)
#   bias 초기화 방식은 zeros_initializer로 변경
cell = tf.contrib.rnn.BasicLSTMCell(num_units=hidden_size, state_is_tuple=True)
outputs, _states = tf.nn.dynamic_rnn(cell, X, dtype=tf.float32)
fc_out=layers.fully_connected(outputs[:,-1],num_classes,activation_fn=tf.nn.relu,
        weights_initializer=layers.xavier_initializer(),biases_initializer=tf.zeros_initializer())
weights=tf.get_variable('weights',shape=[num_classes,out_dim],initializer=layers.xavier_initializer())
bias=tf.get_variable('bias',shape=[out_dim],initializer=tf.zeros_initializer())
logits=tf.matmul(fc_out,weights)+bias
#   prediction값은 softmax를 사용해 학습된 결과를 뽑고
#   cost는 softmax_cross_entropy를 통해 loss값을 구하도록 설정
prediction=tf.nn.softmax(logits,name='prediction')
cost=tf.nn.softmax_cross_entropy_with_logits(logits=logits,labels=Y)
loss=tf.reduce_mean(cost)
loss_summ=tf.summary.scalar("loss",loss)
train = tf.train.AdamOptimizer(learning_rate=learning_rate).minimize(loss)

#   학습된 결과와 Y label을 비교해서 accuracy를 측정   
correct_prediction=tf.equal(tf.argmax(prediction,1),tf.argmax(Y,1))
accuracy=tf.reduce_mean(tf.cast(correct_prediction,tf.float32))

#   텐서보드를 통한 validation의 loss값 확인
merged_summ=tf.summary.merge_all()
writer=tf.summary.FileWriter("./training_logs")
writer.add_graph(sess.graph)

#   학습 내용을 저장(테스트에서 불러오기 위함)
saver=tf.train.Saver()
sess.run(tf.global_variables_initializer())

#   초기 training epoch은 10000 회
#   epoch 당 트레이닝 데이터를 셔플 후 배치로 나눠서 학습
#   epoch 마다 validation loss를 확인하고 저장은 20회마다 하도록 설정
for epoch in range(3177):
    np.random.shuffle(xy_train)
    for i in range(int(16000/batch_size)):
        xy_batch=xy_train[i*batch_size:i*batch_size+batch_size]
        x_batch,y_batch = xy_batch[:,1:-2].reshape(-1,seq_leng,1),xy_batch[:,-2:]

        l,_=sess.run([loss,train],feed_dict={X:x_batch,Y:y_batch})
    print(epoch, l)
    if epoch%20==0:
        saver.save(sess,'./learning_model/training.ckpt')
    
    summary=sess.run(merged_summ,feed_dict={X:x_val,Y:y_val})
    writer.add_summary(summary,global_step=epoch)

Y_real=[]
Y_pred=[]
accu=0
#   테스트도 배치로 나눠서 실행
#   각 테스트 당 실제 데이터와 학습 결과 데이터를 Y_real, Y_pred에 저장
#   0은 snp(1,0) 1은 정상(0,1)을 의미
for i in range(int(np.ceil(6119/batch_size))):
    test_batch=xy_test[i*batch_size:i*batch_size+batch_size]
    test_x,test_y=test_batch[:,1:-2].reshape(-1,seq_leng,1),test_batch[:,-2:]

    accu+=sess.run(accuracy,feed_dict={X:test_x,Y:test_y})

    Y_real+=list(sess.run(tf.argmax(Y,1),feed_dict={Y:test_y}))
    Y_pred+=list(sess.run(tf.argmax(prediction,1),feed_dict={X:test_x}))

print('accuracy : ',accu/int(np.ceil(6119/batch_size)))

#   Y_ret : 3(정상 찾음) 2(정상 못찾음) 1(snp 못찾음) 0(snp 찾음)
#   카운터를 통한 학습 결과를 출력
Y_ret=2*np.array(Y_real)+np.array(Y_pred)
print(Counter(Y_ret))
