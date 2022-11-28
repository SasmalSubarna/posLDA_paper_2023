#!/bin/bash

# path for the python script
f_path=./FES_from_Reweighting.py

#----+++++few notes++++-------#
# derivative option --der not allowed with periodic cv like dihedral
# --stride can't be used with --blocks


# command lines 

# not considering the bias due to walls applied
#python $f_path --colvar ../COLVAR --temp 360.0 --cv ld1 --sigma 0.038 --min -15.0 --max 5.0 --bin 200 --outfile fe_ld1_rew.dat

# fe wrt 6 state lda
#python $f_path --colvar colvar-out.txt --temp 360.0 --cv ld1 --sigma 0.0156 --min -5.0 --max 15.0 --bin 100 --outfile fe_ld1_6_state_rew.dat
python $f_path --colvar colvar-out.txt --temp 360.0 --cv ld1 --sigma 0.0126 --min -5.0 --max 15.0 --bin 100 --outfile fe_ld1_6_state_rew.dat

# considering all the biases applied on the system 
#python $f_path --colvar ../COLVAR --kt 2.577483 --cv cleft_dist,dihedral --sigma 0.03,0.03 --outfile fes_cd+dih_with_walls.dat --skiprows 1000
