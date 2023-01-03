#!/bin/bash

# load all the modules needed
source /scratch/work/hockygroup/software/gromacs-2019.6-plumedSept2020/bin/GMXRC.bash.modules
gmxexe=/scratch/work/hockygroup/software/gromacs-2019.6-plumedSept2020/bin/gmx_mpi

# make topology
mpirun -np 1 $gmxexe pdb2gmx -f hp35_clean.pdb -o hp35_processed.gro -water tip3p -ignh

#==========#
#mpirun -np 1 $gmxexe pdb2gmx -f hp35_clean_no_1LEU.pdb -o hp35_processed_no_1LEU.gro -water tip3p
#mpirun -np 1 $gmxexe pdb2gmx -f hp35_clean_replace.pdb -o hp35_processed_replace.gro -water tip3p 
#===========#

# define a cubic box of length 85A
mpirun -np 1 $gmxexe editconf -f hp35_processed.gro -o hp35_newbox.gro -c -bt cubic -box 8.5 8.5 8.5

# solvate with TIP3P (spc216.gro is thre- point water model and can be used for tip3p)
mpirun -np 1 $gmxexe solvate -cp hp35_newbox.gro -cs spc216.gro -o hp35_solv.gro -p topol.top

# add ions
mpirun -np 1 $gmxexe grompp -f ions.mdp -c hp35_solv.gro -p topol.top -o ions.tpr
mpirun -np 1 $gmxexe genion -s ions.tpr -o hp35_solv_ions.gro -p topol.top -conc 0.04 -pname NA -nname CL -neutral






