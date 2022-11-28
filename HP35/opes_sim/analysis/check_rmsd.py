# import libraries
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as md 

# atoms selection
sel_atoms = np.array([20,22,24,31,33,35,43,45,47,58,60,62,70,72,74,90,92,94,112,114,116,122,124,126,138,140,142,158,160,162,165,167,169,182,184,186,196,198,200,220,222,224,231,233,235,241,243,245,261,263,265,271,273,275,285,287,289,304,306,316,318,320,322,337,339,341,361,363,365,380,382,384,397,399,401,414,416,418,432,434,436,451,453,455,470,472,474,492,494,496,507,509,511,529,531,533,536,538,540,555,557], dtype=int)

sel_atoms -= 1  # mdtraj 0-based indexing 

# reference str
ref_str = md.load("hp35_310K_folded_npt.gro")
print("ref_str:", ref_str)

# load traj 
#traj = md.load("../wt-metad_ld1_hp35_360K.xtc", top=ref_str.topology, stride=10)
traj = md.load("opes_metad_ld1_hp35_360K_bf_8.0_barrier_10.0_wrapped.xtc", top=ref_str.topology, stride=100)
print("traj:", traj)

time = np.arange(traj.n_frames)*5 # in ns  

# superpose
traj_superposed = traj.superpose(reference=ref_str, atom_indices=sel_atoms)

# rmsd 
rmsd = md.rmsd(traj_superposed, ref_str, atom_indices=sel_atoms)*10.0

out = np.array([time,rmsd], dtype=float).T
np.savetxt("rmsd_data_mdtraj.txt", out)

plt.figure(1, figsize=(7,6), dpi=80)
plt.plot(time, rmsd, linestyle="-", lw=2.0, color="navy")
plt.grid(linestyle="--", color="grey", axis="both", which="major")
plt.ylabel("RMSD (angstroms)", fontsize=14)
plt.xlabel("Time (ns)", fontsize=14)
plt.title("OPES_METAD simulation, 360K", fontsize=14)
plt.savefig("rmsd_vs_time.pdf", format="pdf")
plt.show()


