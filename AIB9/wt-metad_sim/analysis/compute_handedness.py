#!/usr/bin/python
import numpy as np
import mdtraj as md
import matplotlib.pyplot as plt

"""
def handedness(trj_file, left_struct='helix_left.xyz', right_struct='helix_right.xyz', top_file='helix.pdb', rmsd_cutoff=0.15):
    trj = md.load(trj_file,top=top_file)
    left_trj = md.load(left_struct,top=top_file)
    right_trj = md.load(right_struct,top=top_file)

    left_rmsd = md.rmsd(trj, left_trj)
    right_rmsd = md.rmsd(trj, right_trj)

    in_left_helix = np.where( (left_rmsd<right_rmsd)*(left_rmsd<rmsd_cutoff))[0]
    in_right_helix = np.where( (left_rmsd>right_rmsd)*(right_rmsd<rmsd_cutoff))[0]

    handedness = np.zeros(len(left_rmsd))
    handedness[in_left_helix] = left_rmsd[in_left_helix]-rmsd_cutoff
    handedness[in_right_helix] = rmsd_cutoff-right_rmsd[in_right_helix]

    return handedness

print("Input pdb handedness: ",handedness('helix.pdb'))
"""

def handedness(r_l, r_r, r_cut):
    """
    r_l = rmsd wrt left
    r_r = rmsd wrt right
    r_cut = cutoff rmsd
    """

    in_l = np.where( (r_l < r_r) * (r_l < r_cut))[0]
    in_r = np.where( (r_l > r_r) * (r_r < r_cut))[0]

    handedness = np.zeros(len(r_l))
    handedness[in_l] = r_l[in_l] - r_cut
    handedness[in_r] = r_cut - r_r[in_r]

    return handedness

# two helix structures
left = md.load("left.gro")
right = md.load("right.gro")
#print(left, right)

# topolgy
topol = left.topology
#print(topol)

# load the traj 
traj = md.load("wt_metad_ld1_aib9_400K_height_0.005_bf_2.0_wrapped.xtc", top=topol, atom_indices=topol.select("protein"), stride=10)
#print(traj.n_frames)

# superposing the trajectory 
traj_l = traj.superpose(left, atom_indices=topol.select("backbone"))
traj_r = traj.superpose(right, atom_indices=topol.select("backbone"))

# rmsd
rmsd_l = md.rmsd(traj_l, left, atom_indices=topol.select("backbone"))*10.0
rmsd_r = md.rmsd(traj_r, right, atom_indices=topol.select("backbone"))*10.0

np.savetxt("rmsd_left.txt", rmsd_l)
np.savetxt("rmsd_right.txt", rmsd_r)

# call handedeness
r_cutoff = 1.50  # in A
handu = handedness(rmsd_l, rmsd_r, r_cutoff) 


# plot
plt.figure(1, figsize=(8,8), dpi=120)
plt.xlabel("RMSD_L (A)", fontsize=16)
plt.ylabel("RMSD_R (A)", fontsize=16)
plt.scatter(rmsd_l, rmsd_r, c=traj.time*1e-3)
cbar = plt.colorbar()
cbar.set_label("Time(ns)", fontsize=16)
plt.savefig("scatter_rmsd_left+right.png")

plt.figure(2, figsize=(8,6), dpi=120)
plt.xlabel("Time (ns)", fontsize=16)
plt.ylabel("RMSD (A)", fontsize=16)
plt.plot(traj.time*1e-3, rmsd_l, label="RMSD_L")
plt.plot(traj.time*1e-3, rmsd_r, label="RMSD_R")
plt.legend()
plt.savefig("rmsds_vs_time.png")

plt.figure(3, figsize=(7,6), dpi=120)
plt.xlabel("Time (ns)", fontsize=16)
plt.ylabel("Handedness", fontsize=16)
plt.title("r_cut=%s A"%str(r_cutoff), fontsize=16)
plt.plot(traj.time*1e-3, handu, color="blue")
plt.savefig("handu_vs_time.png")

plt.show()





























