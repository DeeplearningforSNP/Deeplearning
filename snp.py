import random

input_file = open('hg19_chr21.txt', 'r')

data = input_file
output_file = open('ret.txt', 'w')
snp = open('snp.txt', 'w')
op = ""
agtc = ['a', 'g', 't', 'c', 'A', 'G', 'T', 'C']
data = input_file.read()
rand = []
i = 0

while True:
    while True:
        temp = int(random.random() * 100000000)
        if (temp > 49092490):
            continue

        pos = temp - 7
        a = data[temp].lower()
        if a == 'a' or a == 'g' or a == 'c' or a == 't':
            # print(data[pos])
            if i != 0:
                for j in range(i-1):
                    pos == rand[j]
                    continue
            rand = rand + [pos]
            break
    while True:
        change = random.choice(agtc)
        if (data[temp].lower() != change):
            break

    row = int(rand[i] / 51) + 2
    col = rand[i] % 51 + 1
    wr = "%d(%d/%d) = %c -> %c\n" % ((row - 2) * 50 + col, row, col, data[temp], change)
    output_file.write(wr)
    op = data[:temp] + change + data[temp + 1:]
    i += 1
    if (i == 100):
        break

snp.write(op)

input_file.close()
output_file.close()
snp.close()
