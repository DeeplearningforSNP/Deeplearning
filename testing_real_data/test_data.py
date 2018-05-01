import tensorflow as tf
import numpy as np

batch_size=1000
seq_leng=9

sess=tf.Session()
new_saver=tf.train.import_meta_graph('./learning_model/sum.ckpt.meta')
new_saver.restore(sess,tf.train.latest_checkpoint('./learning_model/'))

graph=tf.get_default_graph()
X=graph.get_tensor_by_name('X:0')
prediction=graph.get_tensor_by_name('prediction:0')

test=open('real_qsum1.csv','r')
max_gene_count=test.read().split('\n')[-1]

ret=open('ret_sum.txt','w')
Y_pred=[]
check_list=[]

test_set=test.read().split('\n')[:-1]
for i in range(len(test_set)):
    data=test_set[i].split(', ')[1:]
    data[8]/=(max_gene_count*2)
    check_list+=[[int(test_est[i].split(', ')[0])]+data[4:8]]
    test_set[i]=data
test_set=np.array(test_set,dtype=np.float32)

for i in range(int(np.ceil(len(test_set)/batch_size))):
    test_batch=test_set[i*batch_size:i*batch_size+batch_size].reshape(-1,seq_leng,1)
    Y_pred+=list(sess.run(tf.argmax(prediction,1),feed_dict={X:test_batch}))

pred_snp=[i for i,x in enumerate(Y_pred) if x==0]
for i in range(len(pred_snp)):
    ret.write("%s"%check_list[i][0])
    if max(check_list[i][1:])==check_list[i][1]:
        ret.write(", A\n")
    elif max(check_list[i][1:])==check_list[i][2]:
        ret.write(", T\n")
    elif max(check_list[i][1:])==check_list[i][3]:
        ret.write(", G\n")
    elif max(check_list[i][1:])==check_list[i][4]:
        ret.write(", C\n")
test.close()
ret.close()
