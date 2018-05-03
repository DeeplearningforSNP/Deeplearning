snp_file=open('sample_result.txt','r')
snp_line=snp_file.read().split('\n')
snp_data=[0]*10000
for i in range(10000):
    snp_data[i]=[int(snp_line[i].split('(')[0])-1]+[snp_line[i].split('>')[1][1]]
learn_file=open('sample_ret.txt','r')
learn_line=learn_file.read().split('\n')
learn_data=[]
for i in range(len(learn_line)-1):
    learn_data+=[learn_line[i].split(', ')]
    learn_data[i][0]=int(learn_data[i][0])
print(len(snp_data),len(learn_data))
snp_dic={'correct':0,'not same':0,'find another':0,'not found':0}
for i in snp_data:
    right=0
    for j in learn_data:
        if i[0]==j[0]:
            right+=1
            if i[1].upper()==j[1]:
                snp_dic['correct']+=1
            else:
                snp_dic['not same']+=1
            break
    if right==0:
        snp_dic['not found']+=1
snp_dic['find another']=len(learn_data)+snp_dic['not found']-len(snp_data)
print(snp_dic)
snp_file.close()
learn_file.close()
