# Import libraries
import numpy as np
import os
import sys
import mdtraj as md
from shapeGMM import _traj_tools as traj_tools
import matplotlib.pyplot as plt
import time 


def kabsch_rotate(mobile, ref):
    correlation = np.dot(mobile.T,ref)
    v, s, w_tr = np.linalg.svd(correlation)
    # enforce determinant of rotation matrix is 1 (not negative 1)
    if np.linalg.det(v) * np.linalg.det(w_tr) < 0.0:
        v[:,-1] *= -1
    rotation = np.dot(v, w_tr)
    #print(rotation)
    # rotate
    return np.dot(mobile,rotation)

def project_trajectory_lda_weighted(trajectory, ref, precision, lda_vector):
    # meta data
    n_frames = trajectory.shape[0]
    n_atoms = trajectory.shape[1]
    # check everything else has matching size
    assert ref.shape[0] == n_atoms, "Atom numbers different in trajectory and reference"
    assert precision.shape[0] == n_atoms, "Atom numbers different in trajectory and precision"
    assert lda_vector.shape[0] == n_atoms*3, "LDA vector dimension do not match trajectory"
    # declare arrays
    projections = np.empty(n_frames,dtype=np.float64)
    # loop over trajectory
    weighted_ref = np.dot(precision, ref)
    for frame in range(n_frames):
        disp = (kabsch_rotate(trajectory[frame],weighted_ref) - ref).flatten()
        projections[frame] = np.dot(disp,lda_vector)
    return projections

def E_lda(R,R0,precision0,v):
    assert R.shape == R0.shape, "Dimensions do not match"
    n_atoms = R.shape[0]
    rot_R = kabsch_rotate(R,np.dot(precision0,R0))
    disp = (rot_R - R0).flatten()
    return np.dot(disp,v.flatten())


def lda_numeric_derivative(R,R0,precision0,v,delta):
    dR = np.empty(R.shape,dtype=np.float64)
    n_atoms = R.shape[0]

    for atom in range(n_atoms):
        for j in range(3):
            
            item = R[atom,j]
            R[atom,j] += delta
            
            dR[atom,j] = E_lda(R,R0,precision0,v)
            
            R[atom,j] = item

    dR -= E_lda(R,R0,precision0,v)
    dR /= delta
    return dR



## read reference and precision
#ref_ = np.loadtxt("hp35_weighted_global_avg.dat")
#prec_ = np.loadtxt("hp35_weighted_global_precision.dat")
#vec_ = np.loadtxt("hp35_weighted_ld1.dat")

ref_ = np.loadtxt("global_avg.txt")
prec_ = np.loadtxt("global_precision.txt")
vec_ = np.loadtxt("lda_coeff_global_6_states.txt")

# path
path_to_files = "./"
#traj_file = path_to_files+"pnas2012-2f4k-360K-protein-000.dcd"
traj_file = path_to_files+"single_frame.dcd"
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


# get lda projection from aligned traj 
lda_projections = project_trajectory_lda_weighted(out, ref_, prec_, vec_)
print(lda_projections)

st = time.time()

numeric = lda_numeric_derivative(out[0,:,:], ref_, prec_, vec_.reshape(101,3),0.00001)

et = time.time()

nt = (et-st)*1e3 # convert to ms 

#print(numeric)


print("nt = ", nt, "ms")

"""
plt.plot(lda_projections)
plt.xlabel("# frame", fontsize=14)
plt.ylabel("LD1", fontsize=14)
plt.grid(axis="both", which="major", linestyle="--")
plt.show()
"""


