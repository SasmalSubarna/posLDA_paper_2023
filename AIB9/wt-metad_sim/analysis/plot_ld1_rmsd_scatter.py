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
plt.xlabel("LD1 (A)", fontsize=16)
plt.ylabel("RMSD_L (A)", fontsize=16)
#plt.scatter(ld1[::100], rmsd_l, c=time[::100], cmap="jet")
plt.scatter(ld1[::100], rmsd_l, c=zeta, cmap="inferno")
cbar = plt.colorbar()
#cbar.set_label("Time (ns)", fontsize=16)
cbar.set_label("$\\zeta$ (radian)", fontsize=16)
plt.savefig("left_rmsd+ld1_scatter.png")


plt.figure(2, figsize=(6,6), dpi=100)
plt.xlabel("LD1 (A)", fontsize=16)
plt.ylabel("RMSD_R (A)", fontsize=16)
#plt.scatter(ld1[::100], rmsd_r, c=time[::100], cmap="jet")
plt.scatter(ld1[::100], rmsd_r, c=zeta, cmap="inferno")
cbar = plt.colorbar()
#cbar.set_label("Time (ns)", fontsize=16)
cbar.set_label("$\\zeta$ (radian)", fontsize=16)
plt.savefig("right_rmsd+ld1_scatter.png")


plt.show()
