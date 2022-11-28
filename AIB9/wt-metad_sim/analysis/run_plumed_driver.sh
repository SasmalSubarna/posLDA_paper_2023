source ./plumed_source_script.sh

# command line
plumed --no-mpi driver --plumed plumed.dat --mf_xtc ../wt-metad_ld1_hp35_360K.xtc
