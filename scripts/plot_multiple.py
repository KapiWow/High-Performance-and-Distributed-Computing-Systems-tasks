#!/usr/bin/python

import os
import sys

exec_name = sys.argv[1]
n_exec = int(sys.argv[2])
line_count = int(sys.argv[3])
count_exec = 10

n_list = [n_exec]

n=n_exec
for i in range(8):
    n_exec = n_exec+n
    n_list.append(n_exec)

res = []

for i in n_list:
    avg = 0
    time_list = [0] * line_count
    print(i)
    for j in range(count_exec): 
        os.system('./' + exec_name + ' ' + str(i))
        file = open("test", "r")
        for k in range(len(time_list)):
            time_list[k] += int(file.readline())
    for j in range(len(time_list)):
        time_list[j] = time_list[j]/count_exec
    res.append(time_list)

print(res)

lines = []
header = "n"
for i in range(len(n_list)):
    header = header + str(i) + " "
lines.append(header)

for i in range(len(n_list)):
    line = str(n_list[i])
    for j in res[i]:
        line += " "
        line += str(j)
    print(line)
    line += '\n'
    lines.append(line)

f = open("data","w+")
f.writelines(lines)

