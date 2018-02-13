import numpy as np
import pymbar

# We need 1) list of restraints, 2) restraint variables


def compute(rest_list, rest_data, method='MBAR'):
    
    #max_n = max([rest_data[i][0].size for i in range(len(rest_data))])
    #k_n = np.zeros([len(rest_data)], np.int32)
    #u_kn = np.zeros([max_n*max_n], np.float64)
    #for pot_win in range(len(rest_data)):
    #    k_n[pot_win] = len(rest_data[pot_win][0])

    #attach_restraints = [rest for rest in rest_list if rest.phase['attach']['force_constants'] is not None]

    kB = 1.381e-23 * 6.022e23 / (4.184 * 1000.0) # Boltzmann constant in kJ/mol/K
    temp = 298.15
    beta = 1/(kB * temp) # beta


    active_rest = [rest.phase['attach']['force_constants'] is not None for rest in rest_list]

    if any(active_rest):
        for i,rest in enumerate(rest_list):
            if active_rest[i]:
                reorder = np.argsort(rest.phase['attach']['force_constants'])
                if rest.continuous_apr:
                    continuous_apr = True
                break
        num_win = len(rest.phase['attach']['force_constants'])

        N_k = np.zeros([num_win], np.int32)
        for i in range(num_win):
            N_k[i] = rest_data['attach'][reorder[i]][0].size
        if continuous_apr:
            N_k[-1] = rest_data['pull'][0][0].size
        N_max = np.max(N_k)

        if method == 'MBAR':
            u_kln = np.zeros([num_win,num_win,N_max], np.float64)
        else:
            u_kn = np.zeros([num_win], np.float64)

        #k = 0
        for pot_idx in range(num_win):
            if method == 'MBAR':
                #if continuous_apr and pot_win == num_win-1:
                #    print(rest_list[0].phase['pull']['force_constants'][0],np.mean(rest_data['pull'][0][0]))
                #else:
                #    print(rest_list[0].phase['attach']['force_constants'][sort_idxs[pot_win]],np.mean(rest_data['attach'][sort_idxs[pot_win]][0]))
                
                pass

#                for crd_idx in range(num_win):
#                    for rest_idx,rest in enumerate(rest_list):
#                        if active_rest[i]:
#                            if rest.mask3 is not None and rest.mask4 is not None:
#                                target = rest.phase['attach']['targets'][reorder[pot_idx]]
#                                bool_list = rest_data['attach'][reorder[crd_idx]][rest_idx] < target - 180.0
#                                rest_data['attach'][reorder[crd_idx]][rest_idx][bool_list] += 360.0
#                                bool_list = rest_data['attach'][reorder[crd_idx]][rest_idx] > target + 180.0
#                                rest_data['attach'][reorder[crd_idx]][rest_idx][bool_list] -= 360.0
#                    u_kln[k,l,0:N_k[k]] += np.sum(
#                        beta*
#                        rest_list[active_rest].phase['attach']['force_constants'][reorder[pot_idx]]*
#                        (rest_data['attach'][reorder[crd_idx]][active_rest] - 
#                        rest_list[active_rest].phase['attach']['targets'][reorder[pot_idx]])**2
#                        )
#        
        
    # If attach exists
        # establish attach_windows based on continuous_apr
        # Figure out which restraints have changing force_constant, store in list? changers
        # Sort attach_windows, force_constant based to be linearly increasing.  Deals with added windows.  Note that added windows will always need to go at the end of a list, out of order.
        # find max data points in rest_data for those windows
        # If full_matrix, u_kn = np.zeros([k*k])
        # else, u_kn = np.zeros([k])
        # winnum = 0
        # for pot_win in range(k):
            
            # if full_matrix:
                # for crd_win in range(k)
                    # u_kn[winnum] = np.sum(beta*fcs*(rest_data[] - targs)**2, axis=1)
            # else:
                #
    

        
# we need to generate the reduced potential energy arrays.
