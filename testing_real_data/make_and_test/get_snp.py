import random

reference = open('hg19_chr21.fa', 'r')
ret_file = open('snp_result.txt', 'w')
snp = open('sample.fa', 'w')
agtc = ['a', 'g', 't', 'c', 'A', 'G', 'T', 'C']
data = reference.read()
rand = []
i = 0

while True:
    while True:
        temp = int(random.random() * 100000000)
        if (temp >= 49092500):
            continue

        a = data[temp].lower()
        if a == 'a' or a == 'g' or a == 'c' or a == 't':
            pos = temp - 7
            d=0
            for j in range(i):
                if pos == rand[j]:
                    d+=1
            if d==0:
                rand = rand + [pos]
                break
    while True:
        change = random.choice(agtc)
        if (data[temp].lower() != change.lower()):
            break

    row = int(rand[i] / 51) + 2
    col = rand[i] % 51 + 1
    wr = "%d(%d/%d) = %c -> %c\n" % ((row - 2) * 50 + col, row, col, data[temp], change)
    ret_file.write(wr)
    data = data[:temp] + change + data[temp + 1:]
    i += 1
    if (i == 10000):
        break

snp.write(data)

reference.close()
ret_file.close()
snp.close()
