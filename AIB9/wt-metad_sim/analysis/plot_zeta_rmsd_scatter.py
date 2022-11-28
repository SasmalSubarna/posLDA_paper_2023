import numpy as np
import matplotlib.pyplot as plt

# load
rmsd_l = np.loadtxt("rmsd_left.txt")
rmsd_r = np.loadtxt("rmsd_right.txt")

zeta = np.load("zeta.npy", allow_pickle=True)
print(zeta.shape)

time, ld1 = np.loadtxt("../COLVAR", usecols=(0,1), unpack=True)
time *= 1e-3 # in ns

plt.figure(1, figsize=(6,6), dpi=100)
plt.xlabel("$\\zeta$ (radian)", fontsize=16)
plt.ylabel("RMSD_L (A)", fontsize=16)
plt.scatter(zeta, rmsd_l, c=ld1[::100], cmap="plasma")
cbar = plt.colorbar()
cbar.set_label("LD1 (A)", fontsize=16)
plt.savefig("zeta+left_rmsd_scatter.png")


plt.figure(2, figsize=(6,6), dpi=100)
plt.xlabel("$\\zeta$ (radian)", fontsize=16)
plt.ylabel("RMSD_R (A)", fontsize=16)
plt.scatter(zeta, rmsd_r, c=ld1[::100], cmap="plasma")
cbar = plt.colorbar()
cbar.set_label("LD1 (A)", fontsize=16)
plt.savefig("zeta+right_rmsd_scatter.png")


plt.show()
