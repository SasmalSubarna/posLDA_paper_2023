import numpy as np
import MDAnalysis as md
from numba import jit
import os

############ subroutines

@jit(nopython=True)
def kabsch_rotate(mobile, ref):
    correlation = np.dot(mobile.T,ref)
    v, s, w_tr = np.linalg.svd(correlation)
    # enforce determinant of rotation matrix is 1 (not negative 1)
    if np.linalg.det(v) * np.linalg.det(w_tr) < 0.0:
        v[:,-1] *= -1
    rotation = np.dot(v, w_tr)
    # rotate
    return np.dot(mobile,rotation)

@jit(nopython=True)
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

### main program

## read global avg and precision and lda vector
ref = np.loadtxt("hp35_weighted_global_avg.dat")
precision = np.loadtxt("hp35_weighted_global_precision.dat")
lda_vector = np.loadtxt("hp35_weighted_ld1.dat")

# create MDAnalysis universe
data_path = '../pnas2012-2f4k-360K-protein/'
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

# compute the LDA projection
lda_projections = project_trajectory_lda_weighted(traj_backbone_101, ref, precision, lda_vector)
# save data
np.savetxt("hp35_ld1_project_global_alignment.dat", lda_projections)

