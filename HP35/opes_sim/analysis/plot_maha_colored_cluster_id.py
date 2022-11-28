import numpy as np
import matplotlib.pyplot as plt

# load
time, d0, d4 = np.loadtxt("maha_data.txt", usecols=(0,1,2), unpack=True)
time *= 50.0*1e-3 # in ns 

# load cluster ids
clusters = np.loadtxt("clusters_predicted.txt")

c4_indx = [i[0] for i in np.argwhere(clusters == 4)]
all_indx_excpt_c4 = [i[0] for i in np.argwhere(clusters != 4)]

d0_for_c4 = d0[c4_indx]
d0_for_all_excpt_c4 = d0[all_indx_excpt_c4]

id_= np.argwhere(d0 == max(d0_for_c4))[0][0]*50.0
print(id_, max(d0))

plt.figure(1, figsize=(7,6), dpi=80)
plt.xlabel("Time (ns)", fontsize=14)
plt.ylabel("Mahalanobis Distance\nw.r.t cluster=0", fontsize=14)
plt.plot(time[all_indx_excpt_c4], d0_for_all_excpt_c4, '.', markersize=1.0)
plt.plot(time[c4_indx], d0_for_c4, '.', ms=1.0, label="c_id = 4")
plt.plot(id_*1e-3, 231.84717099, marker="x", ms=12.0, label="MAX", color="k")
plt.legend(fontsize="large", markerscale=1.0, loc="upper right")
plt.savefig("choose_init_pull.pdf", format="pdf")
plt.show()
