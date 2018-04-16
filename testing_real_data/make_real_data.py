import pysam
import random
import numpy as np

ref_file = open('hg19_chr21.fa', 'r')
ref_data = ref_file.read()
samfile = pysam.AlignmentFile('read.bam', 'rb')
data_file1 = open('real_qsum1.csv', 'w')
data_file2 = open('real_qsum2.csv', 'w')
data_file3 = open('real_qsum3.csv', 'w')
data_file4 = open('real_qsum4.csv', 'w')
data_file5 = open('real_qsum5.csv', 'w')

qual_cost = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6,
             '(': 7, ')': 8, '*': 9, '+': 10, ',': 11, '-': 12, '.': 13,
             '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
             '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27,
             '=': 28, '>': 29, '?': 30, '@': 31, 'A': 32, 'B': 33, 'C': 34,
             'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40, 'J': 41}
gene_cost = {'A': 0, 'T': 1, 'G': 2, 'C': 3}
gene_one_hot = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
result = []
max_total = 0

for pileupcolumn in samfile.pileup('chr21'):
    if pileupcolumn.n==0:
        continue
    if pileupcolumn.pos < 10000000:
        data_file=data_file1
    elif pileupcolumn.pos < 20000000:
        data_file=data_file2
    elif pileupcolumn.pos < 30000000:
        data_file=data_file3
    elif pileupcolumn.pos < 40000000:
        data_file=data_file4
    else:
        data_file=data_file5

    data_file.write("%s, "%pileupcolumn.pos)
    
    pos_in_file = int(pileupcolumn.pos / 50) * 51 + pileupcolumn.pos % 50 + 7
    ref_gene = gene_one_hot[gene_cost[ref_data[pos_in_file].upper()]]
    data_file.write("%s, %s, %s, %s, "%(ref_gene[0],ret_gene[1],ref_gene[2],ref_gene[3]))
    
    gene_count=0
    qual_value = [0, 0, 0, 0]
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            gene = pileupread.alignment.query_sequence[pileupread.query_position]
            qual = pileupread.alignment.qual[pileupread.query_position]
            gene_count+=1
            qual_value[gene_cost[gene.upper()]] += qual_cost[qual]
    for i in range(len(qual_value)):
        qual_value[i]/=(41*gene_count)
    data_file.write("%s, %s, %s, %s, %s\n"%(qual_value[0],qual_value[1],qual_value[2],qual_value[3],gene_count))
    
    max_total = max(max_total, gene_count)
data_file5.write("%s"%max_total)

ret_file.close()
samfile.close()
data_file1.close()
data_file2.close()
data_file3.close()
data_file4.close()
data_file5.close()
