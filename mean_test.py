import pysam
import random

ref_file = open('hg19_chr21.fa', 'r')
ref_data = ref_file.read()

snp_file = open('snp_data.txt', 'r')
snp_data = snp_file.read()
snp_line = snp_data.split('\n')
snp_pos=[0]*10000
samfile = pysam.AlignmentFile('read.bam', 'rb')
train_file = open('real_test_qmean.csv', 'w')

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

for i in range(10000):
    snp_pos[i]=snp_line[i].split('(')[0]

for pileupcolumn in samfile.pileup('chr21'):
    

    gene_count = [0, 0, 0, 0]
    qual_value = [0, 0, 0, 0]
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            gene = pileupread.alignment.query_sequence[pileupread.query_position]
            qual = pileupread.alignment.qual[pileupread.query_position]
            gene_count[gene_cost[gene.upper()]] += 1
            qual_value[gene_cost[gene.upper()]] += qual_cost[qual]

    total_count = sum(gene_count)
    if total_count == 0:
        continue
    max_total = max(max_total, total_count)

    pos_in_file = int(pileupcolumn.pos / 50) * 51 + pileupcolumn.pos % 50 + 7
    ref_gene = ref_data[pos_in_file]
    data = ref_gene
    for num in range(4):
        data += "/" + str(gene_count[num]) + "/" + str(qual_value[num])

    snp = 0
    if pileupcolumn.pos == int(snp_pos[i]) - 1:
        snp = 1

    data += "/" + str(total_count) + "/" + str(snp)
    ref_one_hot = gene_one_hot[gene_cost[ref_gene.upper()]]
    write = "%s, %s, %s, %s, %s, " % (pileupcolumn.pos,ref_one_hot[0],ref_one_hot[1],ref_one_hot[2],ref_one_hot[3])
    data = data.split('/')
    for i in range(1,9,2):
        mean_normal=41
        if int(data[i])!=0:
            mean_normal*=int(data[i])
        write+=("%s, %s, "%(int(data[i])/int(data[9]),int(data[i+1])/mean_normal))
    write += ("%s, %s, %s\n"%(int(data[9])/(max_total*2),data[10],int(not int(data[10]))))
    train_file.write(write)
    
    # pos_dic[data] = pileupcolumn.pos
