import numpy as np
import matplotlib.pyplot as plt

time, ld1 =  np.loadtxt("../COLVAR", usecols=(0,1), unpack=True)
zeta = np.load("zeta.npy")

time *= 1e-3

print(ld1.shape)
print(zeta.shape)

#plt.figure(1)
#plt.plot(time[::100], zeta)


plt.figure(2, figsize=(8,8), dpi=120)
plt.xlabel("LD1", fontsize=16)
plt.ylabel("$\\zeta$", fontsize=16)
plt.scatter(ld1[::100], zeta, c=time[::100], marker='o')
cbar = plt.colorbar()
cbar.set_label("Time (ns)", fontsize=16)
plt.savefig("corrl_ld1_zeta.png")
plt.show()

