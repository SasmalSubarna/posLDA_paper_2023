source ./plumed_source_script.sh

# command line
#plumed --no-mpi driver --plumed plumed.dat --mf_dcd /scratch/work/hockygroup/gmh4/projects/gmm_clustering/data/DESRES-Trajectory_pnas2012-2f4k-360K-protein/pnas2012-2f4k-360K-protein/pnas2012-2f4k-360K-protein-000.dcd

plumed --no-mpi driver --plumed plumed_squared_test.dat --mf_dcd pnas2012-2f4k-360K-protein-000.dcd

#plumed --no-mpi driver --plumed plumed.dat --mf_dcd single_frame.dcd
