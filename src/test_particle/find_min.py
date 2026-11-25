#!/usr/bin/env python3

import os
import sys
import numpy as np
from sys import argv

def main():
    
    eta = float(argv[1])

    middles = []
    AB = []

    filename = str(sys.argv[2])

    f = open(filename, "r")

    for i, line in enumerate(f):
        
        if line == '\n' or line == ' \n':
            continue

        else:
            line = line.split()
            a = float(line[0])
            b = float(line[1])
            
            if line[2]=='nan' or line[3] == 'nan':
                continue
            else:
                mu = float(line[2])
                xt = float(line[3])
                
                AB.append([a, b])
                middles.append(0.5*abs(mu)+0.5*abs(xt))
                #middles.append(1.-np.exp(-mu**2-xt**2))

    count = 0
    while count < 10:
        

        mini = np.min(middles)
        mini_index = middles.index(mini)

        #print("{:.4f} {:.4f} {:.4f} {:.4f} {:.6f}".format(eta, AB[mini_index][0], AB[mini_index][1], 8.*AB[mini_index][0] + 2.*AB[mini_index][1]-9., mini))
        print("{}. A = {:.4f}, B = {:.4f}, PY-distance = {:.4f}, min = {:.6f}".format(count+1,AB[mini_index][0], AB[mini_index][1], 8.*AB[mini_index][0] + 2.*AB[mini_index][1]-9.,mini))
        
        middles.remove(mini)
        
        count += 1

if __name__=='__main__':

    main()
