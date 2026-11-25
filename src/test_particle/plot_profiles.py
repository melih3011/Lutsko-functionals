#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import scienceplots
from sys import argv
from mpl_toolkits.axes_grid1.inset_locator import inset_axes


if __name__=='__main__':
  
  data_mc = np.loadtxt(argv[1])
  data_rf = np.loadtxt(argv[2])
  data_wb = np.loadtxt(argv[3])
  data_wb2 = np.loadtxt(argv[4])

  Z_mc = data_mc[:, 0]
  Z = data_rf[:, 0]

  skip = 2

  plt.style.use('science')
  fig, ax = plt.subplots(figsize=(8, 6))
  
  plt.xticks(fontsize=14)
  plt.yticks(fontsize=14)


  ax.set_xlabel(r'$r/\sigma$', size=18)  
  ax.set_ylabel(r'$g(r)$', size=18)  

  ax.set_xlim(0.9, 3)
  ax.set_ylim(0, 4)
  
  ax.plot(Z_mc[::skip], data_mc[::skip, 1], marker='o', markerfacecolor='none', linestyle='none', label='MC', color='purple')
  ax.plot(Z, data_rf[:, 1], label='RF', color='black')
  ax.plot(Z, data_wb[:, 1], label='WB', color='blue')
  ax.plot(Z, data_wb2[:, 1], label='WBII', color='green')

  ax_inset = inset_axes(ax, width=2.3, height=2.3, bbox_to_anchor=(2.20, 3.75), bbox_transform=ax.transData)
  ax_inset.set_xlim(0.95, 1.15)
  ax_inset.set_ylim(2.2, 3.9)

  ax_inset.plot(Z_mc[::skip], data_mc[::skip], marker='o', markerfacecolor='none', linestyle='none', color='purple')
  ax_inset.plot(Z, data_rf[:, 1], color='black')
  ax_inset.plot(Z, data_wb[:, 1], color='blue')
  ax_inset.plot(Z, data_wb2[:, 1], color='green')
  
  ax.legend(fontsize=14)
  
  plt.savefig("rhos.pdf", dpi=300)
  plt.show()



