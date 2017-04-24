import re
import numpy as np
import subprocess as sp

def write_tleapin(filename='tleap.in', tleaplines=[], unitname='model', pdbfile=None,
                  pbctype=1, bufferval=12.0, waterbox='TIP3PBOX', neutralize=1,
                  countercation='Na+', counteranion='Cl-', addion_residues=None,
                  addion_values=None, saveprefix=None):
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
            for i,res in addion_residues:
                f.write("addionsrand {} {} {}\n".format(unitname,res,addion_values[i]))
        if saveprefix:
            f.write("savepdb {} {}.pdb\n".format(unitname,saveprefix))
            f.write("saveamberparm {} {}.prmtop {}.rst7\n".format(unitname,saveprefix,saveprefix))
        f.write("desc {}\n".format(unitname))
        f.write("quit\n")
    
def countwaters(filename='tleap.in'):
    tleapoutput = sp.Popen('tleap -s -f '+filename, stdout=sp.PIPE, stderr=sp.PIPE, shell=True).stdout.read().splitlines()
    wateradded = 0
    for line in tleapoutput:
        if re.search("^R<WAT ",line):
            wateradded += 1
    return wateradded

def solvate(tleapfile,pdbfile=None, pbctype=1, bufferwater='12.0A', waterbox='TIP3PBOX', neutralize=1,
            countercation='Na+', counteranion='Cl-', addions=None):
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
    bufferval.append(0.0)
    bufferval.append(8.0)
    buffer_iterexp = 0
    wat_added = []

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
                      counteranion=None, addion_residues=None, addion_values=None, saveprefix=None)
        bufferwater = countwaters(filename='tleap_apr_solvate.in')

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

    buffer_predict=0
    if buffer_predict == 1:
        # Do bufferval prediction as conceived by Germano        
        write_tleapin(filename='tleap_apr_solvate.in', tleaplines=tleaplines, unitname=unitname, pdbfile=pdbfile,
                      pbctype=pbctype, bufferval=bufferval[-1], waterbox=waterbox, neutralize=neutralize, countercation=countercation,
                      counteranion=counteranion, addion_residues=addion_residues, addion_values=addion_values, saveprefix=None)
        wat_added.append(countwaters(filename='tleap_apr_solvate.in'))
        prev_bufferval = bufferval
        bufferval += 1.0
        write_tleapin(filename='tleap_apr_solvate.in', tleaplines=tleaplines, unitname=unitname, pdbfile=pdbfile,
                      pbctype=pbctype, bufferval=bufferval, waterbox=waterbox, neutralize=neutralize, countercation=countercation,
                      counteranion=counteranion, addion_residues=addion_residues, addion_values=addion_values, saveprefix=None)
        wat_added.append(countwaters(filename='tleap_apr_solvate.in'))
        cur_bufferval = bufferval
        slope = ((wat_added[-1] - wat_added[-2])/(cur_bufferval - prev_bufferval))
        inter = wat_added[-2] - (slope*prev_bufferval)
        bufferval = (bufferwater - inter)/slope
        if abs(bufferval - prev_bufferval) > 20:
            bufferval = prev_bufferval + 20
        write_tleapin(filename='tleap_apr_solvate.in', tleaplines=tleaplines, unitname=unitname, pdbfile=pdbfile,
                      pbctype=pbctype, bufferval=bufferval, waterbox=waterbox, neutralize=neutralize, countercation=countercation,
                      counteranion=counteranion, addion_residues=addion_residues, addion_values=addion_values, saveprefix=None)
        wat_added.append(countwaters(filename='tleap_apr_solvate.in'))

    print bufferval
    cycle = 0
    buffer_iter = 0
    while cycle < 50:
        write_tleapin(filename='tleap_apr_solvate.in', tleaplines=tleaplines, unitname=unitname, pdbfile=pdbfile,
                      pbctype=pbctype, bufferval=bufferval, waterbox=waterbox, neutralize=neutralize, countercation=countercation,
                      counteranion=counteranion, addion_residues=addion_residues, addion_values=addion_values, saveprefix=None) 
        wat_added.append(countwaters(filename='tleap_apr_solvate.in'))
        print cycle,buffer_iter,":",bufferwater,':',bufferval,buffer_iterexp,wat_added[-2],wat_added[-1]
        cycle += 1
        buffer_iter += 1
        if 0 <= (wat_added[-1] - bufferwater) < 15 or buffer_iterexp < -3:
            print 'Done!'
            break
        elif wat_added[-2] > bufferwater and wat_added[-1] > bufferwater:
            bufferval += -1*(10**buffer_iterexp)
        elif wat_added[-2] > bufferwater and wat_added[-1] < bufferwater:
            if buffer_iter > 1:
                buffer_iterexp -= 1
                buffer_iter = 0
                bufferval += 5*(10**buffer_iterexp)
            else:
                bufferval += 1*(10**buffer_iterexp)
        elif wat_added[-2] < bufferwater and wat_added[-1] > bufferwater:
            if buffer_iter > 1:
                buffer_iterexp -= 1
                buffer_iter = 0
                bufferval += -5*(10**buffer_iterexp)
            else:
                bufferval += -1*(10**buffer_iterexp)
        elif wat_added[-2] < bufferwater and wat_added[-1] < bufferwater:
            bufferval += 1*(10**buffer_iterexp)
        else:
            raise Exception("The bufferval search died due to an unanticipated set of variable values")
      

    if cycle >= 50:
        print "Error message!"
    
    

    # do final run and check

    # print statistics

    # return info for bufferwater, etc to keep them the same.

    #return 'solvated.prmtop'

solvate('tleap.in',bufferwater=15001,pbctype=2)
