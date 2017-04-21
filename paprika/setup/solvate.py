import numpy as np

#def bufferoptimize(bufferwater,wateradded,prev_wateradded,cycles,bufferval,exponent):
#    if (exponent < -3 and wateradded > bufferwater) or (wateradded - bufferwater < 10

def solvate(tleapfile,pdbfile=None, pbctype=1, bufferwater='12.0A', waterbox='TIP3PBOX', neutralize=1,
            countercation='Na+', counteranion='Cl-', addions=None, hmr=0):
    """
    This routine solvates a solute system with a specified amount of water/ions.

    ARGUMENTS

    tleapfile : a fully functioning tleap file which is capable of preparing the system in gas phase. It should load all parameter files that will be necessary for solvation, including the water model and any custom residues.

    pdbfile : if present, the function will search for any loadpdb commands in the tleapfile and replace whatever is there with pdbfile.  This would be used for cases where we have modified PDBs.

    pbctype : the type of periodic boundary conditions. 0 = cubic, 1 = rectangular, 2 = truncated octahedron.

    waterbox : the water box name to use with the solvatebox/solvateoct command. 

    neutralize : 0 = do not neutralize the system, 1 = neutralize the system. the counterions to be used are specified below with 'countercation' and 'counteranion'.

    countercation : a mask to specify neutralizing cations

    counteranion : a mask to specify neturalizing anions

    addions : a list of residue masks and values which indicate how much additional ions to add. The format for the values is as following: if the value is an integer, then add that exact integer amount of ions; if the value is followed by an 'M', then add that amount in molarity;  if 'm', add by molality.  example: ['Mg+',5, 'Cl-',10, 'K+','0.050M']

    """

    unitname = 'model'

    # Read tleapfile
    tleaplines = []
    with open(tleapfile, 'r') as f:
        for line in f.readlines():
            if re.search('loadpdb',line) and pdbfile == None:
                words = line.rstrip().replace('=',' ').split()
                pdbfile = words[2]
                unitname = words[0]
            tleaplines.append(line)
   
    # If bufferwater ends with 'A', meaning it is a buffer distance, we need to get this value
    if str(bufferwater).endswith('A') or :
        with open('tleap_apr_solvate.in', 'w') as f:
        # write a new tleap file with solvation and check results, set integer bufferwater value

    # Process Addions List
    addions = ['Mg+',5, 'Cl-','10m', 'K+','0.050M']
    if addions:
        addion_residues =  []
        addion_values = []
        if len(addions) % 2 == 1:
            raise Exception("The 'addions' list requires an even number of elements. "\
                            "Make sure there is a residue mask followed by a value for "\
                            "each ion to be added")
        for i,txt in enumerate(addions):
            if i % 2 == 0:
                addion_residues.append(txt)
            else:
                if str(txt).endswith('m'):
                    # Figure out number of ions for desired molality
                    addion_values.append(float(txt[:-1])*)
                elif str(txt).endswith('M'):
                    # Figure out number of ions for desired molarity
                else:
                    # must be integer, add to list
        
        print addion_residues
        print addion_values

    # write new tleap file with solvate and addion steps

    # start 8.0, go 9.0.  calc the slope and estimate the right value

    # start optimization near correct value

    # do final run and check

    # print statistics

    # return info for bufferwater, etc to keep them the same.

    return 'solvated.prmtop'

print solvate('tleap.in')
