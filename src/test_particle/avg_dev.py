#!/usr/bin/env python3

import os
import sys
import numpy as np

directory_name = os.getcwd()+'/{}'.format(sys.argv[1])
directory = os.fsencode(directory_name)

def main():
    
    As = []
    Bs = []
    N = len(os.listdir(directory))
    #p = 1./6.
    #q = (1.-2.*p)/(N-2)
    p, q = 1./N, 1./N
    size = 0
    fp = open(directory_name+'/'+os.listdir(directory_name)[0], "r")
    for line in fp:
        if line == ' \n':
            As.append('')
            Bs.append('')
        else:
            line = line.split()
            As.append(float(line[0]))
            Bs.append(float(line[1]))
        size += 1

    print(size)

    mus = np.zeros(size)
    xts = np.zeros(size)

    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        f = open(directory_name+'/'+filename, "r")
        for i, line in enumerate(f):            
            if line == ' \n':
                mus[i] = 0.
                xts[i] = 0.
            else:
                line = line.split()
                if filename.count('0.35')>0 or filename.count('0.40')>0:
                    mus[i] += abs(float(line[2]))*p
                    xts[i] += abs(float(line[3]))*p
                else:
                    mus[i] += abs(float(line[2]))*q
                    xts[i] += abs(float(line[3]))*q

    data = open("averages", "w")

    for i in range(size):
       if As[i] == '':
        data.write("\n")
       else:
        data.write("{:.4f}  {:.4f}  {:.14e} {:.14e}\n".format(As[i], Bs[i], mus[i], xts[i]))

if __name__=='__main__':
    
    main()
