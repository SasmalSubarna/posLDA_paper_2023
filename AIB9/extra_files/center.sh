#!/bin/bash 

# loading all the required modules and stuff
source /scratch/work/hockygroup/software/gromacs-2020.4-plumed2020/bin/GMXRC.bash.modules

# centering the protein in trajectory using PBC condition
mpirun -np 1 /scratch/work/hockygroup/software/gromacs-2020.4-plumed2020/bin/gmx_mpi trjconv -f left.gro -s left.tpr -pbc mol -center -o left_centered.gro

