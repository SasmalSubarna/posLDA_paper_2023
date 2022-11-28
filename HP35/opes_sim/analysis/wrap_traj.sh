#-----------------------------------------------------------------------------#
# This script is used to center the protein molecule after simulations.       #

# Because in gromacs, trajectroy file (.xtc) the protein goes out of the box. #
#-----------------------------------------------------------------------------#


# loading all the required modules and stuff
source /scratch/work/hockygroup/software/gromacs-2020.4-plumed2020/bin/GMXRC.bash.modules

# centering the protein in trajectory using PBC condition
mpirun -np 1 /scratch/work/hockygroup/software/gromacs-2020.4-plumed2020/bin/gmx_mpi trjconv -f ../opes_metad_ld1_hp35_360K_bf_8.0_barrier_10.0.xtc -s md.tpr -pbc mol -center -o centered_traj.xtc 

# eliminating the translational and rotational motion of protein molecule from the centered trajectory
mpirun -np 1 /scratch/work/hockygroup/software/gromacs-2020.4-plumed2020/bin/gmx_mpi trjconv -f centered_traj.xtc -s md.tpr -fit rot+trans -o opes_metad_ld1_hp35_360K_bf_8.0_barrier_10.0_wrapped.xtc

#mpirun -np 1 /scratch/work/hockygroup/software/gromacs-2020.4-plumed2020/bin/gmx_mpi trjconv -f metad_ldaCV+dihedral_actin_100ns.trr -s actin-gatpu_ionized_npt_20ns_every50ps-g5.1.4_5mus.tpr -pbc mol -ur compact -center -dt 60 -o md_60ns_short.xtc

 
