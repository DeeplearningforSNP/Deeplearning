import tensorflow as tf
import numpy as np
from collections import Counter

batch_size=50
seq_leng=9

test=open('real_test.csv','r')
xy_test=test.read().split('\n')
pos_dic={}
test_pos=[]
test_datas=[]
for i in range(len(xy_test)-1):
    pos_dic[int(xy_test[i].split(', ')[0])]=np.array(xy_test[i].split(', ')[1:],dtype=np.float32)
    xy_test[i]=xy_test[i].split(', ')[1:]
xy_test=np.array(xy_test[:-1],dtype=np.float32)
np.random.shuffle(xy_test)

sess=tf.Session()

new_saver=tf.train.import_meta_graph('./learning_model/sum.ckpt.meta')
new_saver.restore(sess,tf.train.latest_checkpoint('./learning_model/'))
#sess.run(tf.global_variables_initializer())

graph=tf.get_default_graph()
X=graph.get_tensor_by_name('X:0')
Y=graph.get_tensor_by_name('Y:0')
prediction=graph.get_tensor_by_name('prediction:0')
accuracy=graph.get_tensor_by_name('accuracy:0')

ret=open('ret_sum.txt','w')
Y_real=[]
Y_pred=[]
accu=0
for i in range(int(np.ceil(4650/batch_size))):
    test_batch=xy_test[i*batch_size:i*batch_size+batch_size]
    test_x,test_y=test_batch[:,:-2].reshape(-1,seq_leng,1),test_batch[:,-2:]

    accu+=sess.run(accuracy,feed_dict={X:test_x,Y:test_y})

    Y_real+=list(sess.run(tf.argmax(Y,1),feed_dict={Y:test_y}))
    Y_pred+=list(sess.run(tf.argmax(prediction,1),feed_dict={X:test_x}))

print('accuracy : ',accu/int(np.ceil(4650/batch_size)))

Y_ret=2*np.array(Y_real)+np.array(Y_pred)
print(Counter(Y_ret))
true_snp=[i for i, x in enumerate(Y_ret) if x==0]
false_snp=[i for i, x in enumerate(Y_ret) if x==1]
#print(len(true_snp),len(false_snp))
c=0
for i in range(len(true_snp)):
    for pos in list(pos_dic.keys()):
        if np.array_equal(xy_test[true_snp[i]],pos_dic[pos]):
            ret.write("%s\n"%pos)
            print(c,i,pos)
#            print(c,pos)
            c+=1
ret.write("\n\n")
for i in range(len(false_snp)):
    for pos in list(pos_dic.keys()):
        if np.array_equal(xy_test[false_snp[i]],pos_dic[pos]):
            ret.write("%s\n"%pos)
test.close()
ret.close()
