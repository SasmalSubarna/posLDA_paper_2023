source ./plumed_source_script.sh

plumed --no-mpi driver --plumed plumed_rog.dat --mf_xtc opes_metad_ld1_hp35_360K_bf_8.0_barrier_10.0_wrapped.xtc

#plumed --no-mpi driver --plumed plumed.dat --mf_xtc opes_metad_ld1_hp35_360K_bf_8.0_barrier_10.0_wrapped.xtc

#plumed --no-mpi driver --plumed plumed_rmsd.dat --mf_xtc opes_metad_ld1_hp35_360K_bf_8.0_barrier_10.0_wrapped.xtc

#plumed --no-mpi driver --plumed plumed_copy.dat --mf_xtc opes_metad_ld1_hp35_360K_bf_8.0_barrier_10.0_wrapped.xtc

#plumed driver --plumed plumed.dat --itrr /scratch/projects/hockygroup/ss12902/new_sims_hlda/metad_hlda/centered_no+rot+trans_metad_hlda_actin_500ns.trr
