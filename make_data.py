import pysam
import random

ref_file=open('hg19_chr21.fa','r')
ref_data=ref_file.read()

snp_file=open('change_data.txt','r')
snp_data=snp_file.read()
snp_line=snp_data.split('\n')
snp_pos=[0]*10000

samfile=pysam.AlignmentFile('gene.bam','rb')
train_file=open('train.csv','w')
validation_file=open('val.csv','w')
test_file=open('test.csv','w')

gene_cost={'A':1,'T':2,'G':3,'C':4}
gene_one_hot=[[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]
cur_num=0
result=[]

max_count=0
for i in range(10000):
    snp_pos[i]=snp_line[i].split('(')[0]
    for pileupcolumn in samfile.pileup('chr21',int(snp_pos[i])-50,int(snp_pos[i])+50):
        if pileupcolumn.pos < int(snp_pos[i])-50 or pileupcolumn.pos > int(snp_pos[i])+50:
            continue
        pos_in_file=int(pileupcolumn.pos/50)*51+pileupcolumn.pos%50+7
        ref_gene=ref_data[pos_in_file]
        gene_count=[0,0,0,0]
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

        total_count/=(max_count*2)
        cur_data=[pileupcolumn.pos]+gene_one_hot[gene_cost[ref_gene.upper()]-1]+gene_count+[total_count,snp]

        dup=0
        if cur_num==0:
            result+=[cur_data]
        else:
            for num in range(len(result)):
                if cur_data[1:]==result[num][1:]:
                    dup+=1
                    break
            if dup==0:
                result+=[cur_data]
        cur_num+=1
random.shuffle(result)
snp_num=0
non_snp=0
for num in range(len(result)):
    if result[num][10]==0:
        if non_snp < 1500:
            train_file.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(
                result[num][0],result[num][1],result[num][2],result[num][3],
                result[num][4],result[num][5],result[num][6],result[num][7],
                result[num][8],result[num][9],result[num][10],int(not result[num][10])))
        elif non_snp < 2810:
            test_file.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(
                result[num][0],result[num][1],result[num][2],result[num][3],
                result[num][4],result[num][5],result[num][6],result[num][7],
                result[num][8],result[num][9],result[num][10],int(not result[num][10])))
        elif non_snp < 3244:
            validation_file.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(
                result[num][0],result[num][1],result[num][2],result[num][3],
                result[num][4],result[num][5],result[num][6],result[num][7],
                result[num][8],result[num][9],result[num][10],int(not result[num][10])))
        non_snp+=1
    else:
        if snp_num < 1500:
            train_file.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(
                result[num][0],result[num][1],result[num][2],result[num][3],
                result[num][4],result[num][5],result[num][6],result[num][7],
                result[num][8],result[num][9],result[num][10],int(not result[num][10])))
        elif snp_num < 1690:
            test_file.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(
                result[num][0],result[num][1],result[num][2],result[num][3],
                result[num][4],result[num][5],result[num][6],result[num][7],
                result[num][8],result[num][9],result[num][10],int(not result[num][10])))
        else:
            validation_file.write("%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n"%(
                result[num][0],result[num][1],result[num][2],result[num][3],
                result[num][4],result[num][5],result[num][6],result[num][7],
                result[num][8],result[num][9],result[num][10],int(not result[num][10])))
        snp_num+=1
    
ref_file.close()
snp_file.close()
samfile.close()
train_file.close()
validation_file.close()
test_file.close()
