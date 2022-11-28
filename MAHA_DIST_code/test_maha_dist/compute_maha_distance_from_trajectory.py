import numpy as np
import MDAnalysis as md
#from numba import jit
import os

############ subroutines

#@jit(nopython=True)
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
    

#@jit(nopython=True)
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

### main program

## read reference and precision
ref = np.loadtxt("hp35_native_cluster_coordinates.dat")
precision = np.loadtxt("hp35_native_cluster_precision.dat")

# create MDAnalysis universe
data_path = '/scratch/work/hockygroup/gmh4/projects/gmm_clustering/data/DESRES-Trajectory_pnas2012-2f4k-360K-protein/pnas2012-2f4k-360K-protein/'
prmtopFileName =  data_path + 'pnas2012-2f4k-360K-protein.pdb'
trajFiles = [data_path + files for files in sorted(os.listdir(data_path)) if files.endswith('.dcd')]
coord = md.Universe(prmtopFileName,trajFiles)
# make atom selection
backbone_selection_101 = "(name C and resid 42) or (name C CA N and not resid 42 76) or (name N and resid 76)"
sel_backbone_101 = coord.select_atoms(backbone_selection_101)
print("Number of atoms in trajectory:", coord.atoms.n_atoms)
print("Number of frames in trajectory:",coord.trajectory.n_frames)
print("Number of atoms being analyzed:",sel_backbone_101.n_atoms)
traj_backbone_101 = np.empty((coord.trajectory.n_frames,sel_backbone_101.n_atoms,3),dtype=float)
count = 0
# read trajectory
for ts in coord.trajectory:
    traj_backbone_101[count,:,:] = sel_backbone_101.positions - sel_backbone_101.center_of_geometry()
    count += 1

# compute the Maha distances
maha_cluster = maha_distance_from_trajectory(traj_backbone_101, ref, precision)
# save data
np.savetxt("hp35_maha_distance_to_native_cluster.dat", maha_cluster)

