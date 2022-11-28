#!/bin/bash

source ./plumed_source_script.sh

# print the value of temp in energy units/ kt
#plumed kt --temp 310 --units kj/mol
plumed kt --temp 400.0 --units kcal/mol
#plumed kt --temp 360.0 --units kcal/mol


# actual one
plumed sum_hills --hills ../HILLS --kt 0.794882 --min -60.0 --max 60.0 --bin 150 --outfile fe_vs_ld1.txt --mintozero



# deposit FE vs. CV after addition every 5000 gaussian hills
#plumed sum_hills --hills ../HILLS --stride 10000 --outfile myfes_ --mintozero --kt 2.577483 
#plumed sum_hills --hills ../HILLS --stride 5000 --outfile myfes_ --mintozero --kt 2.577483


# generating histogram data
#plumed sum_hills --histo ../COLVAR --idw cleft_dist,dihedral --sigma 0.025,0.025 --kt 2.577483 --mintozero

