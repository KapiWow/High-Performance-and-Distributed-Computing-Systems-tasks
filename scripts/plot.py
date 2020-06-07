#!/usr/bin/python

import os
import sys

exec_name = sys.argv[1]
data_name = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])
tick = int(sys.argv[5])
count_exec = int(sys.argv[6])

n_list = range(start, end, tick)


res = []

for i in n_list:
    time = 0
    for j in range(count_exec): 
        os.system('./' + exec_name + ' ' + str(i))
        time += int(file.readline())
    res.append(int(time/count_exec))
print(res)

lines = []
header = "n, time \n"
lines.append(header)

for i in range(len(n_list)):
    line = str(n_list[i])
    line += " " + str(res[i]) + '\n'
    lines.append(line)

f = open("data"+data_name,"w")
f.writelines(lines)
f.close()

os.system("gnuplot " + data_name + ".gnuplot")
