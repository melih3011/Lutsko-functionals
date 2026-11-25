import numpy as np
import os
from sys import argv

eta = float(argv[1])
mix = 0.001
data_name = "data_{:.4f}".format(eta)

os.system("rm {}".format(data_name))

radius = 1.6

#while radius < 20.+0.5:
while radius <= 2.5:

  print("{:.2f}".format(radius))

  #os.system("./spone_wbq3 {:.3f} {:.3f} 0 0 0 {:.6f} >> {}".format(radius, eta, mix, data_name))
  #os.system("./spone_lk {:.3f} {:.3f} 1.3 -1.0 0 {:.6f} >> {}".format(radius, eta, mix, data_name))
  #os.system("./spone_lkwb {:.3f} {:.3f} 1.35 -0.85 0 {:.6f} >> {}".format(radius, eta, mix, data_name))
  os.system("./spone_lkwb2 {:.3f} {:.3f} 2 2 0 {:.6f} >> {}".format(radius, eta, mix, data_name))

  if abs(radius-2.0)<=0.05:
    radius += 0.02
  else: radius += 0.05

