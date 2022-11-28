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
ref_ = np.loadtxt("hp35_native_cluster_coordinates.dat")
prec_ = np.loadtxt("hp35_native_cluster_precision.dat")


# path
path_to_files = "./"
traj_file = path_to_files+"pnas2012-2f4k-360K-protein-000.dcd"
pdb_file = path_to_files+"pnas2012-2f4k-360K-protein.pdb"

#load topology
topol = md.load(pdb_file).topology
#print(topol)

# alpha carbons
ca_atoms = np.array([7,9,10,11,15,16,17,23,24,25,32,33,34,40,41,42,51,52,58,60,61,62,65,66,67,72,73,74,83,84,85,87,88,89,95,96,97,102,103,104,113,114,115,119,120,121,124,125,126,135,136,137,140,141,142,148,149,150,156,157,158,163,164,165,171,172,173,185,186,191,193,194,195,202,203,204,211,212,213,221,222,223,229,230,235,237,238,244,246,247,248,255,256,262,264,265,266,268,269,270,276])

# load trajectory
#traj = md.load(traj_file, top=topol, atom_indices=ca_atoms, stride=100)
traj = md.load(traj_file, top=topol, atom_indices=ca_atoms-1)
#print(traj.n_frames)

data = traj.xyz*10.0
out = np.empty(data.shape)

for i in range(data.shape[0]):
    out[i,:,:] = data[i,:,:] - np.mean(data[i,:,:], axis=0)

maha_dist_data = maha_distance_from_trajectory(out, ref_, prec_)

plt.plot(maha_dist_data)
plt.show()



