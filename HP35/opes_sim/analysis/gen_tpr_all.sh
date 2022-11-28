#!/bin/bash 

source /scratch/work/hockygroup/software/gromacs-2019.6-plumedSept2020/bin/GMXRC.bash.modules.subarna
gmxexe=/scratch/work/hockygroup/software/gromacs-2019.6-plumedSept2020/bin/gmx_mpi

mpirun -np 1 $gmxexe grompp -f md_gen_vel.mdp -c frame_879350.gro -p topol.top -o frame_879350.tpr -maxwarn 1 
