import tensorflow as tf
import numpy as np

batch_size=50
seq_leng=9

sess=tf.Session()
new_saver=tf.train.import_meta_graph('./learning_model/sum.ckpt.meta')
new_saver.restore(sess,tf.train.latest_checkpoint('./learning_model/'))

graph=tf.get_default_graph()
X=graph.get_tensor_by_name('X:0')
prediction=graph.get_tensor_by_name('prediction:0')

test1=open('real_qsum1.csv','r')
test2=open('real_qsum2.csv','r')
test3=open('real_qsum3.csv','r')
test4=open('real_qsum4.csv','r')
max_gene_count=test4.read().split('\n')[-1]

ret=open('ret_sum.txt','w')
Y_pred=[]
pos_dic={}
for i in range(4):
    if i==0:
        testfile=test1
    elif i==1:
        testfile=test2
    elif i==2:
        testfile=test3
    else:
        testfile=test4
    test_set=testfile.read().split('\n')[:-1]
    for j in range(len(test_set)):
        pos=int(test_set[j].split(', ')[0])
        data=test_set[j].split(', ')[1:]+[i]+[j]
        data[8]/=max_gene_count
        pos_dic[pos]=np.array(data,dtype=np.float32)
        test_set[j]=data
    test_set=np.array(test_set,dtype=np.float32)
    np.random_shuffle(test_set)

    for j in range(int(np.ceil(len(test_set)/batch_size))):
        test_batch=test_set[j*batch_size:j*batch_size+batch_size][:-2]
        Y_pred+=list(sess.run(tf.argmax(prediction,1),feed_dict={X:test_batch}))

pred_snp=[i for i,x in enumerate(Y_pred) if x==0]
for i in range(len(pred_snp)):
    for pos in list(pos_dic.keys()):
        if np.array_equal(test_set[pred_snp[i]],pos_dic[pos]):
            ret.write("%s",pos)
            if max(pos_dic[pos][4:8])==pos_dic[pos][4]:
                ret.write(", A\n")
            elif max(pos_dic[pos][4:8])==pos_dic[pos][5]:
                ret.write(", T\n")
            elif max(pos_dic[pos][4:8])==pos_dic[pos][6]:
                ret.write(", G\n")
            elif max(pos_dic[pos][4:8])==pos_dic[pos][7]:
                ret.write(", C\n")
test1.close()
test2.close()
test3.close()
test4.close()
ret.close()
