import numpy as np
import matplotlib.pyplot as plt


# load
time, ld1 = np.loadtxt("../COLVAR", usecols=(0,1), unpack=True)
time *= 1e-3 # in ns

plt.figure(1, figsize=(8,6), dpi=80)
plt.plot(time, ld1, linestyle="-", lw=1.0, color="tab:orange")
plt.grid(linestyle="--", color="grey", axis="both", which="major")
plt.ylabel("LD1 (angstroms)", fontsize=14)
plt.xlabel("Time (ns)", fontsize=14)
plt.title("OPES_METAD simulation, 360K", fontsize=14)
plt.savefig("ld1_vs_time.pdf", format="pdf")
plt.show()

