#!/bin/bash

source /scratch/work/hockygroup/software/gromacs-2019.6-plumedSept2020/bin/GMXRC.bash.modules.subarna
gmxexe=/scratch/work/hockygroup/software/gromacs-2019.6-plumedSept2020/bin/gmx_mpi

# giving wrong structures ! (structures can't be minimized)
#mpirun -np 1 $gmxexe trjconv -f hp35_production_1us_360K_wrapped.xtc -s md.tpr -dt 50000 -sep -vel -o .gro -e 150000
#mpirun -np 1 $gmxexe trjconv -f hp35_production_1us_360K_wrapped.xtc -s md.tpr -dt 50000 -pbc mol -center -ur compact -sep -vel -o .gro -e 150000


# select the centering for entire system always
# notes-  We need to always center the whole system while saving some frames from trajectory. 
# From unwrapped traj use -pbc mol -center -ur compact to place it inside the box perfectly. 
# Always check the coordinates for SOL and PROT to make sure they closely match with benchmark structure or not.
# 

#mpirun -np 1 $gmxexe trjconv -f hp35_production_1us_360K.xtc -s md.tpr -dt 50000 -pbc mol -center -ur compact -sep -vel -o _.gro -e 150000

mpirun -np 1 $gmxexe trjconv -f ../opes_metad_ld1_hp35_360K_bf_8.0_barrier_10.0.xtc -s ../md.tpr -dump 879350 -pbc mol -center -ur compact -o frame_879350.gro
