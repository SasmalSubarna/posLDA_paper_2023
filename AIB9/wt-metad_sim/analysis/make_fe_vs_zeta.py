import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md 


# topology
topol = md.load("left.gro").topology
print(topol)

# load traj
traj = md.load("wt_metad_ld1_aib9_400K_height_0.005_bf_2.0_wrapped.xtc", top=topol, atom_indices=topol.select("protein"))
print(traj)


phi_indices, phi_vals = md.compute_phi(traj)   # it gives only last 8 dihedrals, instead of 9
print(phi_indices)
print(phi_vals.shape)

print(phi_vals[:, 1:6:].shape)
zeta = np.sum(phi_vals[:, 1:6:], axis=1)
np.save("zeta.npy", zeta, allow_pickle=True)

#zeta = np.load("zeta.npy", allow_pickle=True)
#print(zeta.shape)

plt.figure(1)
plt.plot(zeta)
plt.xlabel("# frame", fontsize=16)
plt.ylabel("$\zeta$ (radians)", fontsize=16)
plt.grid(axis="both", which="major", linestyle="--", color="grey")
plt.savefig("zeta_vs_time.png")

nbins=140
bias=  np.loadtxt("../COLVAR", usecols=3)[::10]
print(bias.shape)

plt.figure(2)
hist, bins = np.histogram(zeta, bins=nbins, range=[-7,7], density=True, weights=np.exp(bias/0.794882))
grids = (bins[1:]+ bins[:-1])/2.0
fe = -np.log(hist)
fe -= fe.min()
plt.ylim(0.0,20.0)
plt.xlabel("$\zeta$ (radians)", fontsize=16)
plt.ylabel("FE/ kBT", fontsize=16)
plt.grid(axis="both", which="major", linestyle="--", color="grey")
plt.plot(grids, fe, lw=3.0)
plt.savefig("fe_vs_zeta.png")

plt.show()
