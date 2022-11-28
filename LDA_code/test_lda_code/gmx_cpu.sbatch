#!/bin/bash
#SBATCH --job-name=test_
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --time=10:00
#SBATCH --mem=1GB

# load all the modules needed
#source /scratch/work/hockygroup/software/gromacs-2019.6-plumedSept2020/bin/GMXRC.bash.modules.subarna
#gmxexe=/scratch/work/hockygroup/software/gromacs-2019.6-plumedSept2020/bin/gmx_mpi

source ./plumed_source_script.sh
#srun plumed driver --plumed plumed.dat --mf_dcd pnas2012-2f4k-360K-protein-000.dcd --multi 4
srun plumed driver --plumed plumed.dat --mf_dcd pnas2012-2f4k-360K-protein-000.dcd

# command line
#srun $gmxexe mdrun -s md.tpr -deffnm wt-metad_ld1_hp35_310K -plumed plumed.dat -nsteps 50000000 


