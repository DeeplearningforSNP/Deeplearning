import pysam


ref_file=open('hg19_chr21.fa','r')
ref_data=ref_file.read()

snp_file=open('ret2000.txt','r')
snp_data=snp_file.read()
snp_line=snp_data.split('\n')
snp_pos=[0]*2000

samfile=pysam.AlignmentFile('read.bam','rb')
train_file=open('data_for_train.csv','w')
test_file=open('data_for_test.csv','w')

gene_cost={'A':1,'T':2,'G':3,'C':4}
cur_num=0
result_snp=[]
result_not=[]
max_count=0

def Check(data,result, dup):
    for num in range(len(result)):
        if data[1:] == result[num][1:]:
            dup += 1
            return dup
    return 0

def snpornot(data):
    global result_not
    global result_snp
    if(data[10] == 0):
        result_not += [data]
    else:
        result_snp +=[data]

for i in range(2000):
    snp_pos[i]=snp_line[i].split('(')[0]
    n = int(snp_pos[i])-50
    m = int(snp_pos[i])+50

    for pileupcolumn in samfile.pileup('chr21',n,m):
        gene_count=[0,0,0,0]
        ref_gene =[0,0,0,0]
        for pileupread in pileupcolumn.pileups:
            if type(pileupread.query_position)!=int:
                continue
            gene=pileupread.alignment.query_sequence[pileupread.query_position]
            gene_count[gene_cost[gene.upper()]-1]+=1
        total_count=sum(gene_count)
        max_count=max(max_count,total_count)
        if total_count==0:
            continue
        for num in range(len(gene_count)):
            gene_count[num]/=total_count

        snp=0
        if pileupcolumn.pos==int(snp_pos[i])-1:
            snp=1

        pos_in_file = int(pileupcolumn.pos / 50) * 51 + pileupcolumn.pos % 50 + 7
        ATGC = gene_cost[ref_data[pos_in_file].upper()]
        ref_gene[ATGC-1] = 1
        cur_data=[pileupcolumn.pos]+ref_gene+gene_count+[total_count,snp]
#        print(cur_num,cur_data)

        dup=0
        if cur_num==0:
            snpornot(cur_data)
        else:
            if (snp ==0):
                dup = Check(cur_data, result_not, dup)
                if dup == 0:
                    result_not +=cur_data
            else:
                dup = Check(cur_data, result_snp, dup)
                if dup ==0:
                    result_snp +=cur_data
        cur_num

#        print(result)
#        print(len(result))
#
for num in range(len(result_snp)):
     result_snp[num][9]/=(max_count*2)
     result_not[num][9]/=(max_count*2)
     write = ("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(result_snp[num][0],result_snp[num][1],
             result_snp[num][2],result_snp[num][3],result_snp[num][4],result_snp[num][5],result_snp[num][6],result_snp[num][7],result_snp[num][8],result_snp[num][9], result_snp[num][10]))
     write += ("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(result_not[num][0],result_not[num][1],
             result_not[num][2],result_not[num][3],result_not[num][4],result_not[num][5],result_not[num][6],result_not[num][7],result_not[num][8],result_not[num][9], result_not[num][10]))
     if num < int(len(result_snp)*0.7):
         train_file.write(write)
     else:
         test_file.write(write)

snp_file.close()
samfile.close()
train_file.close()
test_file.close()
