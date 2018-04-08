snp_file=open('snp_data.txt','r')
snp_data=snp_file.read()
snp_line=snp_data.split('\n')
snp_pos=[0]*10000
samtools = open("snp.vcf","r")
haplotype = open("haplotype.vcf","r")
unified = open("unified.vcf","r")
SAM = []
HAP = []
UNI = []

sam = samtools.read()
saml = sam.split("\n")
#32
hap = haplotype.read()
hapl = hap.split("\n")
#29
uni = unified.read()
unil = uni.split("\n")
#33
# for i in range()
# print(saml[28])

# print (unil[len(unil)-2])
for i in range (32, len(saml)-1):
    SAM += [[saml[i].split("\t")[1]]+saml[i].split("\t")[3:5]]

for i in range(29, len(hapl)-1):
    HAP += [[hapl[i].split("\t")[1]]+hapl[i].split("\t")[3:5]]

for i in range(33, len(unil)-1):
    UNI += [[unil[i].split("\t")[1]]+unil[i].split("\t")[3:5]]

# print (SAM[1],HAP[1],UNI[1])
for i in range(10000):
    snp_pos[i]=[snp_line[i].split('(')[0]] + [snp_line[i].split('=')[1][1]]+[snp_line[i].split('>')[1][1]]

# print(snp_pos[2][2])
SAMDic ={'correct':0, 'not same':0, 'find another':0, 'not find':0}
HAPDic ={'correct':0, 'not same':0, 'find another':0, 'not find':0}
UNIDic ={'correct':0, 'not same':0, 'find another':0, 'not find':0}
print (len(snp_pos), len(SAM), len(HAP), len(UNI))
for i in SAM:
    right =0
    for j in range(len(snp_pos)):
        pos = snp_pos[j]
        if i[0] == pos[0]:
            if i[1] == pos[1].upper() and i[2] == pos[2].upper():
                SAMDic['correct'] +=1
                right =1
            else:
                SAMDic['not same'] +=1
                right =1
            break
    if right ==0:
        SAMDic['find another'] +=1
SAMDic['not find'] = len(snp_pos)-len(SAM)+SAMDic['find another']

print ("Samtools : ",SAMDic)
for i in HAP:
    right =0
    for j in range(len(snp_pos)):
        pos = snp_pos[j]
        if i[0] == pos[0]:
            if i[1] == pos[1].upper() and i[2] == pos[2].upper():
                HAPDic['correct'] +=1
                right =1
            else:
                HAPDic['not same'] +=1
                right =1
            break
    if right ==0:
        HAPDic['find another'] +=1
HAPDic['not find'] = len(snp_pos)-len(HAP)+HAPDic['find another']

print("Haplotype : ",HAPDic)
for i in UNI:
    right =0
    for j in range(len(snp_pos)):
        pos = snp_pos[j]
        if i[0] == pos[0]:
            if i[1] == pos[1].upper() and i[2] == pos[2].upper():
                UNIDic['correct'] +=1
                right =1
            else:
                UNIDic['not same'] +=1
                right =1
            break
    if right ==0:
        UNIDic['find another'] +=1
UNIDic['not find'] = len(snp_pos)-len(UNI)+UNIDic['find another']
print("Unified : ",UNIDic)