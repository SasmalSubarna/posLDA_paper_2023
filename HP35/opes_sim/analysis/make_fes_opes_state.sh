#!/bin/bash 

f_name=./FES_from_State.py


#-----++++few notes++++-------#
# --cv is not available here 
# --sigma not needed

# command line
#python $f_name --temp 360.0 --state ../STATES --min -15.0 --max 5.0 --bin 200 --outfile fe_ld1_state.dat --der 
#python $f_name --temp 360.0 --state ../STATES --min -15.0 --max 5.0 --bin 100 --outfile fe_ld1_state.dat --der 

#python $f_name --temp 360.0 --state ../state_2.5us.dat --min -15.0 --max 5.0 --bin 100 --outfile fe_ld1_state.dat --der 
python $f_name --temp 360.0 --state ../state_2.5us.dat --min -20.0 --max 8.0 --bin 150 --outfile fe_ld1_state.dat --der 

#python $f_name --kt 2.577483 --state state_copy.data --all_stored --outfile fe_hlda_state.dat --der 

