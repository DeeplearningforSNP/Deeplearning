import tensorflow as tf
import numpy as np

batch_size=50
seq_leng=9
xy_test=[]
test1=open('real_qsum1.csv','r')
test2=open('real_qsum2.csv','r')
test3=open('real_qsum3.csv','r')
test4=open('real_qsum4.csv','r')
test5=open('real_qsum5.csv','r')
xy_test+=test1.read().split('\n')
xy_test+=test2.read().split('\n')
xy_test+=test3.read().split('\n')
xy_test+=test4.read().split('\n')
xy_test+=test5.read().split('\n')
pos_dic={}
for i in range(len(xy_test)-1):
    pos_dic[int(xy_test[i].split(', ')[0])]=np.array(xy_test[i].split(', ')[1:],dtype=np.float32)
    xy_test[i]=xy_test[i].split(', ')[1:]
xy_test=np.array(xy_test[:-1],dtype=np.float32)
np.random.shuffle(xy_test)

sess=tf.Session()

new_saver=tf.train.import_meta_graph('./learning_model/sum.ckpt.meta')
new_saver.restore(sess,tf.train.latest_checkpoint('./learning_model/'))

graph=tf.get_default_graph()
X=graph.get_tensor_by_name('X:0')
prediction=graph.get_tensor_by_name('prediction:0')

ret=open('ret_sum.txt','w')
Y_pred=[]
accu=0
for i in range(int(np.ceil(4650/batch_size))):
    test_batch=xy_test[i*batch_size:i*batch_size+batch_size]
    Y_pred+=list(sess.run(tf.argmax(prediction,1),feed_dict={X:test_batch}))
pred_snp=[i for i,x in enumerate(Y_pred) if x==0]
for i in range(len(pred_snp)):
    for pos in list(pos_dic.keys()):
        if np.array_equal(xy_test[pred_snp[i]],pos_dic[pos]):
            ret.write("%s",pos)
            if max(pos_dic[pos][5:9])==pos_dic[pos][5]:
                ret.write(", A\n")
            elif max(pos_dic[pos][5:9])==pos_dic[pos][6]:
                ret.write(", T\n")
            elif max(pos_dic[pos][5:9])==pos_dic[pos][7]:
                ret.write(", G\n")
            elif max(pos_dic[pos][5:9])==pos_dic[pos][8]:
                ret.write(", C\n")
test1.close()
test2.close()
test3.close()
test4.close()
test5.close()
ret.close()
