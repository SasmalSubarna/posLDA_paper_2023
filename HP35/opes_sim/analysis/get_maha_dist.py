# Import libraries
import numpy as np
import os
import sys
import mdtraj as md
from shapeGMM import _traj_tools as traj_tools
import matplotlib.pyplot as plt

def maha_dist2_with_align(mobile, ref, precision):
    correlation = np.dot(mobile.T,np.dot(precision,ref))
    v, s, w_tr = np.linalg.svd(correlation)
    # enforce determinant of rotation matrix is 1 (not negative 1)
    if np.linalg.det(v) * np.linalg.det(w_tr) < 0.0:
        v[:,-1] *= -1
    rotation = np.dot(v, w_tr)
    # rotate
    mobile_prime = np.dot(mobile,rotation)
    # compute displacement and return Mahalanobis distance
    disp = mobile_prime - ref
    return np.trace(np.dot(disp.T,np.dot(precision,disp)))


def maha_distance_from_trajectory(positions,ref,precision):
    # meta data
    n_frames = positions.shape[0]
    n_atoms = positions.shape[1]
    assert ref.shape[0] == n_atoms, "Reference has different number of atoms as trajectory"
    assert precision.shape[0] == n_atoms, "Precision has different number of atoms as trajectory"
    # declare arrays
    maha_distances = np.empty(n_frames)
    weighted_ref = np.dot(precision,ref)
    # loop over trajectory
    for frame in range(n_frames):
        maha_distances[frame] = maha_dist2_with_align(positions[frame],ref, precision)
    # return distances
    return np.sqrt(maha_distances)


## read reference and precision
ref_ = np.loadtxt("cluster4_avg.txt")
prec_ = np.loadtxt("cluster4_prec.txt")


# path
path_to_files = "./"
traj_file = path_to_files+"opes_metad_ld1_hp35_360K_bf_8.0_barrier_10.0_wrapped.xtc"
pdb_file = path_to_files+"npt.gro"

#load topology
topol = md.load(pdb_file).topology
#print(topol)

# alpha carbons
ca_atoms = np.array([20,22,24,31,33,35,43,45,47,58,60,62,70,72,74,90,92,94,112,114,116,122,124,126,138,140,142,158,160,162,165,167,169,182,184,186,196,198,200,220,222,224,231,233,235,241,243,245,261,263,265,271,273,275,285,287,289,304,306,316,318,320,322,337,339,341,361,363,365,380,382,384,397,399,401,414,416,418,432,434,436,451,453,455,470,472,474,492,494,496,507,509,511,529,531,533,536,538,540,555,557])


# load trajectory
#traj = md.load(traj_file, top=topol, atom_indices=ca_atoms, stride=100)
traj = md.load(traj_file, top=topol, atom_indices=ca_atoms-1)
#print(traj.n_frames)

data = traj.xyz*10.0
out = np.empty(data.shape)

# align
for i in range(data.shape[0]):
    out[i,:,:] = data[i,:,:] - np.mean(data[i,:,:], axis=0)

maha_dist_data = maha_distance_from_trajectory(out, ref_, prec_)
np.savetxt('test.txt', maha_dist_data)

plt.plot(maha_dist_data)
plt.show()



