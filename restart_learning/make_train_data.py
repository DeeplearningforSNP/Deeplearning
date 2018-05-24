import pysam
import random

ref_file = open('hg19_chr21.fa', 'r')
ref_data = ref_file.read()

snp_file = open('read.txt', 'r')
snp_data = snp_file.read().split('\n')
snp_pos = [0] * 10000

samfile = pysam.AlignmentFile('read.bam', 'rb')

train_file = open('learn_train.csv', 'w')
validation_file = open('learn_val.csv', 'w')
test_file = open('learn_test.csv', 'w')

qual_cost = {'!': 0, '"': 1, '#': 2, '$': 3, '%': 4, '&': 5, "'": 6,
             '(': 7, ')': 8, '*': 9, '+': 10, ',': 11, '-': 12, '.': 13,
             '/': 14, '0': 15, '1': 16, '2': 17, '3': 18, '4': 19, '5': 20,
             '6': 21, '7': 22, '8': 23, '9': 24, ':': 25, ';': 26, '<': 27,
             '=': 28, '>': 29, '?': 30, '@': 31, 'A': 32, 'B': 33, 'C': 34,
             'D': 35, 'E': 36, 'F': 37, 'G': 38, 'H': 39, 'I': 40, 'J': 41}
gene_one_hot = {'A':[1, 0, 0, 0], 'T':[0, 1, 0, 0], 'G':[0, 0, 1, 0], 'C':[0, 0, 0, 1]}

result = []
pos_dic = {}
max_total = 0

for i in range(10000):
    snp_pos[i] = snp_data[i].split('(')[0]
snp_pos.sort()

for i in range(10000):
    for pileupcolumn in samfile.pileup('chr21', int(snp_pos[i]) - 10, int(snp_pos[i]) + 10):
        if (pileupcolumn.pos < int(snp_pos[i]) - 10 or pileupcolumn.pos >= int(snp_pos[i]) + 10) or pileupcolumn.n == 0:
            continue
        here_is_snp = 0
        for j in range(i - 10, i + 10):
            if j>=0 and j<10000 and pileupcolumn.pos!=int(snp_pos[i])-1 and pileupcolumn.pos==int(snp_pos[j])-1:
                here_is_snp = 1
                break
        if here_is_snp == 1:
            continue
        gene_count = 0
        qual_value = {'A':0, 'T':0, 'G':0, 'C':0}
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_del and not pileupread.is_refskip:
                gene = pileupread.alignment.query_sequence[pileupread.query_position]
                qual = pileupread.alignment.qual[pileupread.query_position]
                gene_count += 1
                qual_value[gene.upper()] += qual_cost[qual]
        if gene_count == 0:
            continue
        for k in qual_value.keys():
            qual_value[k]/=(41*gene_count)
        qual_data=[qual_value['A'],qual_value['T'],qual_value['G'],qual_value['C']]

        max_total = max(max_total, gene_count)

        pos_in_file = int(pileupcolumn.pos / 50) * 51 + pileupcolumn.pos % 50 + 7
        ref_gene = ref_data[pos_in_file]
        data = ref_gene
        for num in range(4):
            data += "/" + str(qual_data[num])

        snp = 0
        if pileupcolumn.pos == int(snp_pos[i]) - 1:
            snp = 1

        data += "/" + str(gene_count) + "/" + str(snp)
        pos_dic[data] = pileupcolumn.pos

data_list = list(pos_dic.keys())
random.shuffle(data_list)

def write_data(fname, gene_pos, gene_one_hot, data_set, max_total):
    ref_one_hot = gene_one_hot[data_set[0].upper()]
    fname.write("%s, %s, %s, %s, %s, " % (
        gene_pos, ref_one_hot[0], ref_one_hot[1], ref_one_hot[2], ref_one_hot[3]))
    fname.write("%s, %s, %s, %s, " % (data_set[1], data_set[2], data_set[3], data_set[4]))
    fname.write("%s, %s, %s\n" % (int(data_set[5]) / (max_total*2), data_set[6], int(not int(data_set[6]))))

snp_num = 0
non_snp = 0
for num in range(len(data_list)):
    data_set = data_list[num].split('/')
    gene_pos=pos_dic[data_list[num]]
    if int(data_set[6]) == 0:
        if non_snp < 8000:
            write_data(train_file, gene_pos, gene_one_hot, data_set, max_total)
        elif non_snp < 9600:
            write_data(validation_file, gene_pos, gene_one_hot, data_set, max_total)
        elif non_snp < 14600:
            write_data(test_file, gene_pos, gene_one_hot, data_set, max_total)
        non_snp += 1
    else:
        if snp_num < 8000:
            write_data(train_file, num, gene_one_hot, data_set, max_total)
        elif snp_num < 8400:
            write_data(validation_file, gene_pos, gene_one_hot, data_set, max_total)
        else:
            write_data(test_file, gene_pos, gene_one_hot, data_set, max_total)
        snp_num += 1

ref_file.close()
snp_file.close()
samfile.close()
train_file.close()
validation_file.close()
test_file.close()
