import re
import numpy as np
import subprocess as sp

def write_tleapin(filename='tleap.in', tleaplines=[], unitname='model', pdbfile=None,
                  pbctype=1, bufferval=12.0, waterbox='TIP3PBOX', neutralize=1,
                  countercation='Na+', counteranion='Cl-', addion_residues=None,
                  addion_values=None, removewat=None, saveprefix=None):
    # Add argument info .... ugh duplicates solvate
    # Any suggestions for this function?  so. many. arguments.  Put everything in a dict?
    with open(filename, 'w') as f:
        for line in tleaplines:
            f.write(line)
        if pbctype == 0:
            f.write("solvatebox {} {} {:0.5f} iso\n".format(unitname,waterbox,bufferval))
        elif pbctype == 1:
            f.write("solvatebox {} {} {{10.0 10.0 {:0.5f}}}\n".format(unitname,waterbox,bufferval))
        elif pbctype == 2:
            f.write("solvateoct {} {} {:0.5f} iso\n".format(unitname,waterbox,bufferval))
        else:
            raise Exception("Incorrect pbctype value provided: "+pbctype+". Only 0, 1, and 2 are valid")
        if neutralize == 1:
            f.write("addionsrand {} {} 0\n".format(unitname,countercation))
            f.write("addionsrand {} {} 0\n".format(unitname,counteranion))
        if addion_residues:
            for i,res in enumerate(addion_residues):
                f.write("addionsrand {} {} {}\n".format(unitname,res,addion_values[i]))
        if removewat:
            for watnum in removewat:
                f.write("remove {} {}.{}\n".format(unitname,unitname,watnum))
        if saveprefix:
            f.write("savepdb {} {}.pdb\n".format(unitname,saveprefix))
            f.write("saveamberparm {} {}.prmtop {}.rst7\n".format(unitname,saveprefix,saveprefix))
        f.write("desc {}\n".format(unitname))
        f.write("quit\n")
    
def countresidues(filename='tleap.in', returnlist=None):
    # Add argument info
    tleapoutput = sp.Popen('tleap -s -f '+filename, stdout=sp.PIPE, stderr=sp.PIPE, shell=True).stdout.read().splitlines()
    wateradded = 0
    reslist = []
    for line in tleapoutput:
        r = re.search("^R<WAT (.*)>",line)
        if r:
            wateradded += 1
            if returnlist == 'WAT':
                reslist.append(r.group(1))
    if returnlist == 'ALL':
        reslist = {}
        for line in tleapoutput:
            r = re.search("^R<(.*) (.*)>",line)
            if r:
                if r.group(1) not in reslist:
                    reslist[r.group(1)] = []
                reslist[r.group(1)].append(r.group(2))
        return reslist
    elif returnlist == 'WAT':
        return reslist
    else:
        return wateradded


def solvate(tleapfile, pdbfile=None, pbctype=1, bufferwater='12.0A', waterbox='TIP3PBOX', neutralize=1,
            countercation='Na+', counteranion='Cl-', addions=None, saveprefix='solvated'):
    """
    This routine solvates a solute system with a specified amount of water/ions.

    ARGUMENTS

    tleapfile : a fully functioning tleap file which is capable of preparing the system in gas phase. It should load all parameter files that will be necessary for solvation, including the water model and any custom residues. Assumes the final conformations of the solute are loaded via PDB.

    pdbfile : if present, the function will search for any loadpdb commands in the tleapfile and replace whatever is there with pdbfile.  This would be used for cases where we have modified PDBs.

    pbctype : the type of periodic boundary conditions. 0 = cubic, 1 = rectangular, 2 = truncated octahedron.  If rectangular, only the z-axis buffer can be manipulated; the x- and y-axis will use a 10 Ang buffer.

    waterbox : the water box name to use with the solvatebox/solvateoct command. 

    neutralize : 0 = do not neutralize the system, 1 = neutralize the system. the counterions to be used are specified below with 'countercation' and 'counteranion'.

    countercation : a mask to specify neutralizing cations

    counteranion : a mask to specify neturalizing anions

    addions : a list of residue masks and values which indicate how much additional ions to add. The format for the values is as following: if the value is an integer, then add that exact integer amount of ions; if the value is followed by an 'M', then add that amount in molarity;  if 'm', add by molality.  example: ['Mg+',5, 'Cl-',10, 'K+','0.050M']

    """

    unitname = 'model'
    addion_residues =  []
    addion_values = []
    bufferval = []
    bufferval.append(28.0)
    buffer_iterexp = 1
    wat_added = []
    wat_added.append(0.0)

    # Read tleapfile
    tleaplines = []
    with open(tleapfile, 'r') as f:
        for line in f.readlines():
            if re.search('loadpdb',line) and pdbfile == None:
                words = line.rstrip().replace('=',' ').split()
                if pdbfile == None:
                    pdbfile = words[2]
                unitname = words[0]
                tleaplines.append("{} = loadpdb {}\n".format(unitname,pdbfile))
            if not re.search("^\s*addions|^\s*addions2|^\s*addionsrand|^\s*desc|"\
                             "^\s*quit|^\s*solvate|loadpdb|^\s*save",line,re.IGNORECASE):
                tleaplines.append(line)
   
    # If bufferwater ends with 'A', meaning it is a buffer distance, we need to get this value.
    if str(bufferwater).endswith('A'):
        bufferval.append(float(bufferwater[:-1]))
        write_tleapin(filename='tleap_apr_solvate.in', tleaplines=tleaplines, unitname=unitname, pdbfile=pdbfile,
                      pbctype=pbctype, bufferval=bufferval[-1], waterbox=waterbox, neutralize=0, countercation=None,
                      counteranion=None, addion_residues=None, addion_values=None, removewat=None, saveprefix=None)
        bufferwater = countresidues(filename='tleap_apr_solvate.in')

    # Process Addions List
    #addions = ['Mg+',5, 'Cl-','10m', 'K+','0.050M']
    if addions:
        if len(addions) % 2 == 1:
            raise Exception("Error: The 'addions' list requires an even number of elements. "\
                            "Make sure there is a residue mask followed by a value for "\
                            "each ion to be added")
        for i,txt in enumerate(addions):
            if i % 2 == 0:
                addion_residues.append(txt)
            else:
                if str(txt).endswith('m') or str(txt).endswith('M'):  # Temporary treat m and M identical
                    # Figure out number of ions for desired molality
                    addion_values.append(float(txt[:-1])*float(bufferwater)*0.018)
                else:
                    addion_values.append(txt)

    # First adjust wat_added by changing the bufferval
    cycle = 0
    buffer_iter = 0
    while cycle < 50:
        write_tleapin(filename='tleap_apr_solvate.in', tleaplines=tleaplines, unitname=unitname, pdbfile=pdbfile,
                      pbctype=pbctype, bufferval=bufferval[-1], waterbox=waterbox, neutralize=neutralize, countercation=countercation,
                      counteranion=counteranion, addion_residues=addion_residues, addion_values=addion_values, removewat=None, saveprefix=None) 
        wat_added.append(countresidues(filename='tleap_apr_solvate.in'))
        print cycle,buffer_iter,":",bufferwater,':',bufferval[-1],buffer_iterexp,wat_added[-2],wat_added[-1]
        cycle += 1
        buffer_iter += 1
        if 0 <= (wat_added[-1] - bufferwater) < 12 or buffer_iterexp < -3:
            # Possible failure mode: if the tolerance here is very small (0 < () < 1),
            # the loop can exit with bufferval that adds fewer waters than bufferwater
            print 'Done!'
            break
        ### Possible location of a switch to adjust the bufferval by polynomial fit approach.
        elif wat_added[-2] > bufferwater and wat_added[-1] > bufferwater:
            bufferval.append(bufferval[-1] + -1*(10**buffer_iterexp))
        elif wat_added[-2] > bufferwater and wat_added[-1] < bufferwater:
            if buffer_iter > 1:
                buffer_iterexp -= 1
                buffer_iter = 0
                bufferval.append(bufferval[-1] + 5*(10**buffer_iterexp))
            else:
                bufferval.append(bufferval[-1] + 1*(10**buffer_iterexp))
        elif wat_added[-2] < bufferwater and wat_added[-1] > bufferwater:
            if buffer_iter > 1:
                buffer_iterexp -= 1
                buffer_iter = 0
                bufferval.append(bufferval[-1] + -5*(10**buffer_iterexp))
            else:
                bufferval.append(bufferval[-1] + -1*(10**buffer_iterexp))
        elif wat_added[-2] < bufferwater and wat_added[-1] < bufferwater:
            bufferval.append(bufferval[-1] + 1*(10**buffer_iterexp))
        else:
            raise Exception("The bufferval search died due to an unanticipated set of variable values")
      

    if cycle >= 50:
        raise Exception("Automatic adjustment of the buffer value was unable to converge on a solution with sufficient tolerance")
    elif wat_added[-1] - bufferwater < 0:
        raise Exception("Automatic adjustment of the buffer value resulted in fewer waters added than targeted by bufferwater. Try increasing the tolerance in the above loop")
    else:
        watover = 0
        cycle = 0
        while wat_added[-1] != bufferwater or cycle == 0:
            watover += wat_added[-1] - bufferwater ### Note I don't think there should be water removal errors, but if so, this loop and '+=' method is an attempt to fix.
            watlist = countresidues(filename='tleap_apr_solvate.in', returnlist='WAT')
            if watover == 0:
                removewat = None
            else:
                removewat = watlist[-1*watover:]
            write_tleapin(filename='tleap_apr_solvate.in', tleaplines=tleaplines, unitname=unitname, pdbfile=pdbfile,
                      pbctype=pbctype, bufferval=bufferval[-1], waterbox=waterbox, neutralize=neutralize, countercation=countercation,
                      counteranion=counteranion, addion_residues=addion_residues, addion_values=addion_values, removewat=removewat, saveprefix=saveprefix)
            reslist = countresidues(filename='tleap_apr_solvate.in', returnlist='ALL')
            wat_added.append(len(reslist['WAT']))
            for key,value in sorted(reslist.items()):
                print key,len(value)
            cycle += 1
            if cycle >= 10:
                raise Exception("Solvation failed due to an unanticipated problem with water removal")


solvate('tleap.in',bufferwater=2003,pbctype=1,addions = ['K+',5,'BR',5])
