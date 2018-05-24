import random

ref_file=open('hg19_chr21.fa','r')
data=ref_file.read()

snp_file=open('read.fa','w')
snp_info=open('read.txt','w')

change={'A':'tTgGcC', 'T':'aAgGcC', 'G':'aAtTcC', 'C':'aAtTgG'}
num=0
gene_pos=[]

while True:
    temp=int(random.random()*100000000)
    if temp >= 49092500:
        continue

    gene=data[temp].upper()
    if not gene in 'ATGC':
        continue

    pos=temp-7
    dup=0
    for i in range(num):
        if pos==gene_pos[i]:
            dup+=1
    if dup>0:
        continue
            
    gene_pos+=[pos]
    snp=random.choice(change[gene])

    row=int(gene_pos[num]/51)+2
    col=gene_pos[num]%51+1
    snp_info.write("%d(%d/%d) = %c -> %c\n"%((row-2)*50+col,row,col,data[temp],snp))
    data=data[:temp]+snp+data[temp+1:]
    num+=1
    if num==10000:
        break
snp_file.write("%s"%data)

ref_file.close()
snp_file.close()
snp_info.close()
