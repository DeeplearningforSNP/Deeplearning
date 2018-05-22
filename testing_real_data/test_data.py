import tensorflow as tf
import numpy as np

batch_size=5000
seq_leng=9

sess=tf.Session()
new_saver=tf.train.import_meta_graph('./learning_model/training.ckpt.meta')
new_saver.restore(sess,tf.train.latest_checkpoint('./learning_model/'))

graph=tf.get_default_graph()
X=graph.get_tensor_by_name('X:0')
keep_prob=graph.get_tensor_by_name('keep_prob:0')
prediction=graph.get_tensor_by_name('prediction:0')

test=open('real_qsum1.csv','r')
line = open('real_qsum4.csv','r')
test_data=test.read()
test_set=test_data.split('\n')[:-1]
max_gene_count=line.read().split('\n')[-1]

ret=open('ret_test1.txt','w')
Y_pred=[]
check_list=[]

for i in range(len(test_set)):
    data=test_set[i].split(', ')[1:]
    data[8]=int(data[8])/(int(max_gene_count)*2)
    check_list+=[[int(test_set[i].split(', ')[0])]+data[4:8]]
    test_set[i]=data
test_set=np.array(test_set,dtype=np.float32)

for i in range(int(np.ceil(len(test_set)/batch_size))):
    test_batch=test_set[i*batch_size:i*batch_size+batch_size].reshape(-1,seq_leng,1)
    Y_pred+=list(sess.run(tf.argmax(prediction,1),feed_dict={X:test_batch,keep_prob:1}))

pred_snp=[i for i,x in enumerate(Y_pred) if x==0]
for i in range(len(pred_snp)):
    ret.write("%s"%check_list[pred_snp[i]][0])
    if max(check_list[pred_snp[i]][1:])==check_list[pred_snp[i]][1]:
        ret.write(", A\n")
    elif max(check_list[pred_snp[i]][1:])==check_list[pred_snp[i]][2]:
        ret.write(", T\n")
    elif max(check_list[pred_snp[i]][1:])==check_list[pred_snp[i]][3]:
        ret.write(", G\n")
    elif max(check_list[pred_snp[i]][1:])==check_list[pred_snp[i]][4]:
        ret.write(", C\n")
test.close()
ret.close()
