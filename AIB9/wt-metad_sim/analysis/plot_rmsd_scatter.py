import numpy as np
import matplotlib.pyplot as plt


# load 
rmsd_l = np.loadtxt("rmsd_left.txt")
rmsd_r = np.loadtxt("rmsd_right.txt")

time, ld1 = np.loadtxt("../COLVAR", usecols=(0,1), unpack=True)
time *= 1e-3

zeta = np.load("zeta.npy", allow_pickle=True)

plt.figure(1, figsize=(6,6), dpi=100)
plt.xlabel("RMSD_L (A)", fontsize=16)
plt.ylabel("RMSD_R (A)", fontsize=16)
plt.scatter(rmsd_l, rmsd_r, c=zeta, cmap="viridis")
plt.plot(rmsd_l, rmsd_l, '--', color="k", lw=2.0)
cbar = plt.colorbar()
cbar.set_label("$\\zeta$ (radian)", fontsize=16)
plt.savefig("rmsds_scatter_color_zeta.png")


plt.figure(2, figsize=(6,6), dpi=100)
plt.xlabel("RMSD_L (A)", fontsize=16)
plt.ylabel("RMSD_R (A)", fontsize=16)
plt.scatter(rmsd_l, rmsd_r, c=ld1[::100], cmap="viridis")
plt.plot(rmsd_l, rmsd_l, '--', color="k", lw=2.0)
cbar = plt.colorbar()
cbar.set_label("LD1 (A)", fontsize=16)
plt.savefig("rmsds_scatter_color_ld1.png")

plt.show()
