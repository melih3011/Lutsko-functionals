#!/usr/bin/env python3

import numpy as np
import scipy.integrate as integrate
import matplotlib.pyplot as plt
from sys import argv

pi = np.pi

def main():
  
  data = np.loadtxt(str(argv[1]))
  
  R = 100
  N = len(data)
  
  rhob = 3.*0.35/(4.*pi)
  
  rvals = 2.*data[:, 0]+1.
  rho = rhob*data[:, 1]
  
  rhor2 = np.array([(rho[i+R]-rhob)*rvals[i]**2 for i in range(N-R)])
  
  print(1.+4.*pi*(integrate.simpson(rhor2, rvals[:-R])-8./3.*rhob))
  

if __name__ == "__main__":
  
  main()
