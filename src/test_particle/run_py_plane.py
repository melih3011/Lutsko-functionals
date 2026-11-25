#!/usr/bin/env python3

import numpy as np
import os


if __name__=='__main__':

  deta = 0.01
  a = 0.1
  etas = np.arange(0.01,0.50,deta)
  
  os.system("make clean && make lk")

  for eta in etas:

    print("{:.3f}".format(eta))
    
    if (eta < 0.35): a = 0.1
    elif (eta < 0.40): a = 0.05
    elif (eta < 0.45): a = 0.01

    os.system("./spone_lk 1.0 {:.3f} 0.0 0.0 2 {:.3f} > temp".format(eta, a, eta))
    os.system("./find_min.py {:.3f} temp >> mins.dat".format(eta))
