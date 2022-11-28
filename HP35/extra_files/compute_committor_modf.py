def b_before_a(label_trajectory, state_a=0, state_b=1):
    # save teh committor labels
    import numpy as np
    locations_a = np.where(label_trajectory==state_a)[0]
    locations_b = np.where(label_trajectory==state_b)[0]
    #print(len(locations_a), locations_a[:10:])
    #print(len(locations_b), locations_b[:10:])
    #final_known_state = max(locations_a.max(),locations_b.max())
    #print(final_known_state)
    if len(locations_a) == 0:
        #final_known_state = locations_b.max()
        locations_a= np.array([-1])
    if len(locations_b) == 0:
        #final_known_state = locations_a.max()
        locations_b = np.array([-1])
    #else:
    #    final_known_state = max(locations_a.max(),locations_b.max())
    
    final_known_state = max(locations_a.max(), locations_b.max())

    locations_a_iter = iter(locations_a)
    locations_b_iter = iter(locations_b)

    first_a = next(locations_a_iter)
    first_b = next(locations_b_iter)

    #print(first_a)
    #print(first_b)

    commitor = np.zeros( len(label_trajectory) )
    #c = 0
    current_index = 0
    while first_a > -1 and first_b > -1:
        #c += 1
        if first_b < first_a: 
            commitor[current_index:first_b+1] = 1
            current_index = first_b
            first_b = next(locations_b_iter, -1)
        else:
            current_index = first_a + 1
            first_a = next(locations_a_iter, -2)
        #if c<200:
        #    print("first_a,first_b=", first_a, first_b)

    if first_a < 0:
        commitor[current_index:max(locations_b)+1] = 1

    commitor[final_known_state+1:] = -1

    return commitor

"""
if __name__ == "__main__":
    import mdtraj as md
    import numpy as np
    import sys
    from matplotlib import pylab as plt
    folding_trajectory = md.load("helix_folding_eps6.0.dcd",top="helix_template.pdb")[:-61]
    print(folding_trajectory.n_frames)
    helix_left = md.load('helix_left.xyz',top="helix_template.pdb")
    helix_right = md.load('helix_right.xyz',top="helix_template.pdb")
    
    folding_trj_comp_left = folding_trajectory.superpose(helix_left)
    folding_trj_comp_right = folding_trajectory.superpose(helix_right)
    rmsd_left = md.rmsd(folding_trj_comp_left,helix_left)
    rmsd_right = md.rmsd(folding_trj_comp_right,helix_right)

    trj_labels = 2*np.ones(len(folding_trajectory))
    trj_labels[rmsd_left<0.04] = 0
    trj_labels[rmsd_right<0.04] = 1
    
    # save the traj labels
    np.savetxt("trj_labels.txt", trj_labels, fmt="%1.1f")

    commitor = b_before_a(trj_labels)

    # save the committor labels
    np.savetxt("committor_labels.txt", commitor, fmt="%1.1f")

    #np.savetxt(sys.stdout, np.append(trj_labels.reshape(-1,1),commitor.reshape(-1,1), axis=-1) ,fmt='%d')

"""
