import numpy as np

def solvate(tleapfile, pbctype=1, bufferwater='12.0A', waterbox='TIP3PBOX', neutralize=1,
            countercation='Na+', counteranion='Cl-', addions=None, hmr=0):
    """
    This routine solvates a solute system with a specified amount of water/ions.

    ARGUMENTS

    tleapfile : a fully functioning tleap file which is capable of preparing the system in gas phase. It should load all parameter files that will be necessary for solvation, including the water model and any custom residues.

    pbctype : the type of periodic boundary conditions. 0 = cubic, 1 = rectangular, 2 = truncated octahedron.

    waterbox : the water box name to use with the solvatebox/solvateoct command. 

    neutralize : 0 = do not neutralize the system, 1 = neutralize the system. the counterions to be used are specified below with 'countercation' and 'counteranion'.

    countercation : a mask to specify neutralizing cations

    counteranion : a mask to specify neturalizing anions

    addions : a list of residue masks and values which indicate how much additional ions to add. The format for the values is as following: if the value is an integer, then add that exact integer amount of ions; if the value is followed by an 'M', then add that amount in molarity;  if 'm', add by molality.  example: ['Mg+',5, 'Cl-',10, 'K+','0.050M']

    """

    return 'solvated.prmtop'

print solvate('tleap.in')
