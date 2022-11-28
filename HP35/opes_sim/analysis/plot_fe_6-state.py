import numpy as np
import matplotlib.pyplot as plt

# load
grid, fe = np.loadtxt("fe_ld1_6_state_rew.dat", usecols=(0,1), unpack=True)
fe /= 4.184 # in kcal/mol unit

# load DESRES data (used same params, bins=100, range=[-5,15])
grid_ref, fe_ref = np.loadtxt("fe_vs_ld1_6_states_global_alignment_bins_100_range-5,15.txt", usecols=(0,1), unpack=True)

plt.figure(1, figsize=(7,6), dpi=80)
plt.title("6-State LD1: obtained from 2-state LD1", fontsize=18)
plt.xlabel("LD 1 (angstroms)", fontsize=16)
plt.ylabel("FE (kcal/mol)", fontsize=16)
plt.ylim(0.0,12.0)
plt.grid(linestyle="--", color="grey", axis="both", which="major")
plt.plot(grid_ref, fe_ref, label="ref. FE", lw=3.0, color="k")
plt.plot(grid, fe, label="FE", lw=3.0, color="tab:red")
plt.legend(ncol=1, loc="upper right", fontsize="x-large")
plt.savefig("fe_vs_ld1_6-state.pdf", format="pdf")
plt.show()

