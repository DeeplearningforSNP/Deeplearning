import tensorflow as tf
import pysam
import numpy as np

seq_leng=9
batch_size=100000

ref_file = open('hg19_chr21.fa', 'r')
ref_data = ref_file.read()
samfile = pysam.AlignmentFile('snp.bam', 'rb')

check_list=[]
data_list=[]

sess=tf.Session()
new_saver=tf.train.import_meta_graph('./learning_model/training.ckpt.meta')
new_saver.restore(sess,tf.train.latest_checkpoint('./learning_model/'))

graph=tf.get_default_graph()
X=graph.get_tensor_by_name('X:0')
prediction=graph.get_tensor_by_name('prediction:0')

ret=open('ret_test.txt','w')
Y_pred=[]
max_count=0

qual_cost = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6,
             '(': 7, ')': 8, '*': 9, '+': 10, ',': 11, '-': 12, '.': 13,
             '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
             '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27,
             '=': 28, '>': 29, '?': 30, '@': 31, 'A': 32, 'B': 33, 'C': 34,
             'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40, 'J': 41}
gene_one_hot = {'A':[1, 0, 0, 0], 'T':[0, 1, 0, 0], 'G':[0, 0, 1, 0], 'C':[0, 0, 0, 1]}

for pileupcolumn in samfile.pileup('chr21'):
    if pileupcolumn.n==0:
        continue
    gene_count=0
    qual_value = {'A':0, 'T':0, 'G':0, 'C':0}
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            gene = pileupread.alignment.query_sequence[pileupread.query_position]
            qual = pileupread.alignment.qual[pileupread.query_position]
            gene_count+=1
            qual_value[gene.upper()] += qual_cost[qual]
    if gene_count==0:
        continue
    max_count=max(max_count,gene_count)
    for i in qual_value.keys():
        qual_value[i]/=(41*gene_count)
    qual_data=[qual_value['A'],qual_value['T'],qual_value['G'],qual_value['C']]
   
    pos_in_file = int(pileupcolumn.pos / 50) * 51 + pileupcolumn.pos % 50 + 7
    ref_gene = gene_one_hot[ref_data[pos_in_file].upper()]
   
    check_list+=[pileupcolumn.pos,qual_data]
    data_list+=[ref_gene+qual_data+[gene_count]]

for i in range(int(np.ceil(len(data_list)/batch_size))):
    for j in range(i*batch_size,i*batch_size+batch_size):
        data_list[j][8]/=(max_count*2)
    batch_list=np.array(data_list[i*batch_size,i*batch_size+batch_size],
            dtype=np.float32).reshape(-1,seq_leng,1)
    Y_pred+=list(sess.run(rf.argmax(prediction,1),feed_doct={X:batch_list}))

pred_snp=[i for i,x in enumerate(Y_pred) if x==0]
for i in range(len(pred_snp)):
    ret.write("%s"%check_list[pred_snp[i]][0])
    if max(check_list[pred_snp[i]][1])=='A':
        ret.write(", A\n")
    elif max(check_list[pred_snp[i]][1])=='T':
        ret.write(", T\n")
    elif max(check_list[pred_snp[i]][1])=='G':
        ret.write(", G\n")
    elif max(check_list[pred_snp[i]][1])=='C':
        ret.write(", C\n")

ref_file.close()
samfile.close()
ret.close()
