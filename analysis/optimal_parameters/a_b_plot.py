import numpy as np
import matplotlib.pyplot as plt
import scienceplots

# -------------------
# Provided data (example placeholders)
# Replace these with your actual data
# -------------------

eta = np.array([np.pi/6.*0.2, np.pi/6.*0.4, 0.314159, 0.340339, 0.366519, 0.392699, 0.418879, 0.445059])

#A_MSE = np.array([1.49734, 1.43365, 1.43088, 1.43567, 1.4357, 1.43856, 1.43905, 1.43856])
#A_BV  = np.array([1.47746, 1.38782, 1.35043, 1.33493, 1.32637, 1.31801, 1.32152, 1.30641])

#B_MSE = np.array([-1.51952, -1.61297, -1.7261, -1.83507, -1.82512, -1.91735, -1.89518, -1.97157])
#B_BV  = np.array([-1.31726, -1.33656, -1.30923, -1.32785, -1.28509, -1.32971, -1.33039, -1.34318])

#A_MSE = np.array([1.53556, 1.54949, 1.63309, 1.65804, 1.68061, 1.70353, 1.72725, 1.74639])
#A_BV = np.array([1.51365, 1.45965, 1.46169, 1.45596, 1.45797, 1.46044, 1.47642, 1.47194])

#B_MSE = np.array([-1.38614, -1.52845, -1.77155, -1.92771, -1.95642, -2.09966, -2.11881,-2.24855])
#B_BV = np.array([-0.749039, -0.835034, -0.850284, -0.890756, -0.852691, -0.923799, \
#-0.935022, -0.973172])

A_MSE = np.array([1.56892, 1.57056, 1.63683, 1.65791, 1.67838, 1.69847, 1.72107, 1.73776])
A_BV  = np.array([1.54855, 1.52578, 1.55817, 1.55907, 1.56797, 1.57709, 1.60055, 1.60148])

B_MSE = np.array([-1.55622, -1.73252, -1.92625, -2.07385, -2.08377, -2.21776, -2.2168, \
-2.33626])
B_BV  = np.array([-1.3018, -1.45141, -1.49417, -1.54798, -1.51, -1.5955, -1.61006, \
-1.65833])


# -------------------
# Plot
# -------------------

fig, axes = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
plt.style.use('science')

# ---- Left panel: A ----
axes[0].plot(
    eta, A_MSE,
    color='blue', marker='o', linestyle='-',
    label=r'$A$ MSE'
)
axes[0].plot(
    eta, A_BV,
    color='red', marker='s', linestyle='--',
    label=r'$A$ BV'
)

axes[0].set_xlabel(r'$\eta$')
axes[0].legend()
axes[0].grid(False)

# ---- Right panel: B ----
axes[1].plot(
    eta, B_MSE,
    color='blue', marker='o', linestyle='-',
    label=r'$B$ MSE'
)
axes[1].plot(
    eta, B_BV,
    color='red', marker='s', linestyle='--',
    label=r'$B$ BV'
)

axes[1].set_xlabel(r'$\eta$')
axes[1].legend()
axes[1].grid(False)

# -------------------
# Layout
# -------------------
plt.tight_layout()
plt.savefig("c2_lkwb2.pdf", dpi=300)
plt.show()

