snp_file=open('snp_data','r')
snp_line=snp_file.read().split('\n')
snp_data=[0]*10000
for i in range(10000):
    snp_data[i]=[snp_line[i].split('(')[0]]+[snp_line[i].split('>')[1][1]]
learn_file=open('ret_sum.txt','r')
learn_line=learn_file.read().split('\n')
learn_data=[]
for i in range(len(learn_data)-1):
    learn_data+=[learn_data[i].split(', ')]
snp_dic={'correct':0,'not same':0,'find anothe':0,'not find':0}
for i in snp_pos:
    right=0
    for j in learn_pos:
        if i[0]==j[0]:
            right+=1
            if i[1].upper()==j[1]:
                snp_dic['correct']+=1
            else:
                snp_dic['not same']+=1
            break
    if right==0:
        snp_dic['find another']+=1
snp_dic['not find']=len(snp_data)-len(learn_data)+snp_dic['find another']
print(snp_dic)
