import numpy as np
import matplotlib.pyplot as plt


# load
time, ld1 = np.loadtxt("ld-data.txt", usecols=(0,1), unpack=True)
time *= 50*1e-3 # in ns

plt.figure(1, figsize=(8,6), dpi=80)
plt.plot(time, ld1, linestyle="-", lw=1.0, color="tab:green")
plt.grid(linestyle="--", color="grey", axis="both", which="major")
plt.ylabel("LD1 (angstroms)", fontsize=14)
plt.xlabel("Time (ns)", fontsize=14)
plt.title("OPES_METAD simulation, 360K", fontsize=14)
plt.savefig("ld1_vs_time_6-state.pdf", format="pdf")
plt.show()

