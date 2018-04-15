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
    gene_count=0
    qual_value = [0, 0, 0, 0]
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            gene = pileupread.alignment.query_sequence[pileupread.query_position]
            qual = pileupread.alignment.qual[pileupread.query_position]
            gene_count+=1
            qual_value[gene_cost[gene.upper()]] += qual_cost[qual]

    if gene_count == 0:
        continue
    for i in range(len(qual_value)):
        qual_value[i]/=(41*gene_count)
    max_total = max(max_total, gene_count)

    pos_in_file = int(pileupcolumn.pos / 50) * 51 + pileupcolumn.pos % 50 + 7
    ref_gene = ref_data[pos_in_file]
    result+=[[pileupcolumn.pos]+ref_gene+qual_value+[gene_count]]

re_leng=len(result)
re_part=int(np.ceil(re_leng/5))
for i in range(re_leng):
    result[i][9]/=(max_total*2)
    if i<re_part:
        data_file1.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(
        result[i][0],result[i][1],result[i][2],result[i][3],result[i][4],
        result[i][5],result[i][6],result[i][7],result[i][8],result[i][9]))
    elif i<re_part*2:
        data_file2.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(
        result[i][0],result[i][1],result[i][2],result[i][3],result[i][4],
        result[i][5],result[i][6],result[i][7],result[i][8],result[i][9]))
    elif i<re_part*3:
        data_file3.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(
        result[i][0],result[i][1],result[i][2],result[i][3],result[i][4],
        result[i][5],result[i][6],result[i][7],result[i][8],result[i][9]))
    elif i<re_part*4:
        data_file4.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(
        result[i][0],result[i][1],result[i][2],result[i][3],result[i][4],
        result[i][5],result[i][6],result[i][7],result[i][8],result[i][9]))
    else:
        data_file5.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(
        result[i][0],result[i][1],result[i][2],result[i][3],result[i][4],
        result[i][5],result[i][6],result[i][7],result[i][8],result[i][9]))
ret_file.close()
samfile.close()
data_file1.close()
data_file2.close()
data_file3.close()
data_file4.close()
data_file5.close()
