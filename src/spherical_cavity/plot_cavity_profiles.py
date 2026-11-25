

import numpy as np
import matplotlib.pyplot as plt
import scienceplots


plt.rcParams['font.family'] = 'Arial'

rf = np.loadtxt('sp_lk_4.00_0.32_1.30_-1.00')
wb = np.loadtxt('sp_lkwb_4.00_0.32_1.35_-0.85')
wb2 = np.loadtxt('sp_lkwb_II_4.00_0.32_1.25_-0.20')


fig = plt.figure(figsize=(8, 6))

plt.style.use('science')
plt.xlabel(r'$r/\sigma$', size=19, fontweight='bold')
plt.ylabel(r'$\rho(r)\sigma^3$', size=19, fontweight='bold')
plt.xticks(size=16)
plt.yticks(size=16)

plt.xlim(0, 1.55)
plt.ylim(0, 0.47)

plt.plot(rf[:, 0], rf[:, 1], linewidth=2, color=(0, 0.8, 0.3), label='LK(1.3, -1.00)')
plt.plot(wb[:, 0], wb[:, 1], linewidth=2, color='red', linestyle='--', label='LK-WB(1.35, -0.85)')
plt.plot(wb2[:, 0], wb2[:, 1], linewidth=2, color='blue', linestyle=':', label='LK-WBII(1.25, -0.20)')

plt.grid(0)
legend = plt.legend(fontsize=18)

legend.get_texts()[0].set_color((0, 0.8, 0.3))
legend.get_texts()[1].set_color('red')
legend.get_texts()[2].set_color('blue')

plt.savefig('output.pdf', dpi=300)
plt.show()



