import numpy as np
import pymbar

# We need 1) list of restraints, 2) restraint variables


def compute_reduced_pot(rest_list, rest_data):
    
    max_n = max([rest_data[i][0].size for i in range(len(rest_data))])
    u_kn = np.zeros([max_n*max_n], np.float64)
    for pot_win in range(len(rest_data)):
        for crd_win in range(len(rest_data)):
            
        
# we need to generate the reduced potential energy arrays.
