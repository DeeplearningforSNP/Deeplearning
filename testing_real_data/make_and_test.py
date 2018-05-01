import tensorflow as tf
import pysam
import random
import numpy as np

seq_leng=10
batch_size=1000

ref_file = open('hg19_chr21.fa', 'r')
ref_data = ref_file.read()
samfile = pysam.AlignmentFile('read.bam', 'rb')

data_list=[]

sess=tf.Session()
new_saver=tf.train.import_meta_graph('./learning_model/sum.ckpt.meta')
new_saver.restore(sess,tf.train.latest_checkpoint('./learning_model/'))

graph=tf.get_dafault_graph()
X=graph.get_tensor_by_name('X:0')
prediction=graph.get_tensor_by_nemr('prediction:0')

ret=open('ret_test.txt','w')
Y_pred=[]
pos_list=[]

qual_cost = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6,
             '(': 7, ')': 8, '*': 9, '+': 10, ',': 11, '-': 12, '.': 13,
             '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
             '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27,
             '=': 28, '>': 29, '?': 30, '@': 31, 'A': 32, 'B': 33, 'C': 34,
             'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40, 'J': 41}
gene_cost = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
gene_one_hot = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]

for pileupcolumn in samfile.pileup('chr21'):
    if pileupcolumn.n==0:
        continue
    gene_count=0
    qual_value = [0, 0, 0, 0]
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            gene = pileupread.alignment.query_sequence[pileupread.query_position]
            qual = pileupread.alignment.qual[pileupread.query_position]
            gene_count+=1
            qual_value[gene_cost[gene.upper()]] += qual_cost[qual]
    if gene_count==0:
        continue
    for i in range(len(qual_value)):
        qual_value[i]/=(41*gene_count)
   
    pos_in_file = int(pileupcolumn.pos / 50) * 51 + pileupcolumn.pos % 50 + 7
    ref_gene = gene_one_hot[gene_cost[ref_data[pos_in_file].upper()]]
   
    gene_power=0
    while True:
        if gene_count<1:
            break
        gene_count/=10
        gene_power+=.1

    check_list+=[[pileupcolumn.pos]+qual_value]
    data_list+=[ref_gene+qual_value+[gene_count,gene_power]]
    if len(data_list)==batch_size:
        data_list=np.array(data_list,dtype=np.float32).reshape(-1,seq_leng,1)
        Y_pred+=list(sess.run(tf.argmax(prediction,1),feed_dict={X:data_list}))
        data_list=[]

if len(data_list)!=0:
    data_list=np.array(data_list,dtype=np.float32)
    Y_pred+=list(sess.run(tf.argmax(prediction,1),feed_dict={X:data_list}))

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

ret_file.close()
samfile.close()
ret.close()
