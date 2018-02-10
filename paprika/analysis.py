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
    attach_bool = [rest.phase['attach']['force_constants'] is not None for rest in rest_list]

    if any(attach_bool):
        for i,rest in enumerate(rest_list):
            if attach_bool[i]:
                sort_idxs = np.argsort(rest.phase['attach']['force_constants'])
                if rest.continuous_apr:
                    continuous_apr = True
                break
        NumWin = len(sort_idxs)
        #if continuous_apr:
        #    NumWin += 1

        if method == 'MBAR':
            u_kn = np.zeros([NumWin*NumWin], np.float64)
        else:
            u_kn = np.zeros([NumWin], np.float64)

        k = 0
        print(NumWin,sort_idxs)
        for pot_win in range(NumWin):
            if method == 'MBAR':
                print(pot_win, sort_idxs[pot_win])
                if continuous_apr and pot_win == NumWin-1:
                    print(rest_list[0].phase['pull']['force_constants'][0],np.mean(rest_data['pull'][0][0]))
                else:
                    print(rest_list[0].phase['attach']['force_constants'][sort_idxs[pot_win]],np.mean(rest_data['attach'][sort_idxs[pot_win]]))
                #for crd_win in range(NumWin):
                #    u_kn[k] =
        
        
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
