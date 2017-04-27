import re as re
import os as os
import numpy as np
import subprocess as sp
import logging as log


def write_tleapin(
        filename='tleap.in',
        directory=None,
        tleaplines=[],
        unit='model',
        pbctype=1,
        bufferval=12.0,
        waterbox='TIP3BOX',
        neutralize=1,
        counter_cation='Na+',
        counter_anion='Cl-',
        addion_residues=None,
        addion_values=None,
        removewat=None,
        saveprefix=None):
    """
    Write a `tleap` input file.
    """

    with open(directory + '/' + filename, 'w') as f:
        for line in tleaplines:
            f.write(line)
        if pbctype == 0:
            f.write("solvatebox {} {} {:0.5f} iso\n".format(
                unit, waterbox, bufferval))
        elif pbctype == 1:
            f.write("solvatebox {} {} {{10.0 10.0 {:0.5f}}}\n".format(
                unit, waterbox, bufferval))
        elif pbctype == 2:
            f.write("solvateoct {} {} {:0.5f} iso\n".format(
                unit, waterbox, bufferval))
        else:
            raise Exception("Incorrect pbctype value provided: " +
                            pbctype + ". Only 0, 1, and 2 are valid")
        if neutralize == 1:
            f.write("addionsrand {} {} 0\n".format(unit, counter_cation))
            f.write("addionsrand {} {} 0\n".format(unit, counter_anion))
        if addion_residues:
            for i, res in enumerate(addion_residues):
                f.write("addionsrand {} {} {}\n".format(
                    unit, res, addion_values[i]))
        if removewat:
            for watnum in removewat:
                f.write("remove {} {}.{}\n".format(unit, unit, watnum))
        if saveprefix:
            f.write("savepdb {} {}.pdb\n".format(unit, saveprefix))
            f.write("saveamberparm {} {}.prmtop {}.rst7\n".format(
                unit, saveprefix, saveprefix))
        f.write("desc {}\n".format(unit))
        f.write("quit\n")


def countresidues(filename='tleap.in', directory='.', returnlist=None):
    """
    Run and parse `tleap` output and return the number of residues added.
    returnlist can be ALL for all residuse or WAT for waters.
    """
    # This could probably be replaced with Popen.communicate(), but not sure.
    tleapoutput = sp.Popen('tleap -s -f ' + filename, stdout=sp.PIPE,
                           stderr=sp.PIPE, shell=True, cwd=directory).stdout.read().splitlines()
    wateradded = 0
    reslist = []

    for line in tleapoutput:
        line = line.decode('utf-8')
        r = re.search("^R<WAT (.*)>", line)
        if r:
            wateradded += 1
            if returnlist == 'WAT':
                reslist.append(r.group(1))
    if returnlist == 'ALL':
        reslist = {}
        for line in tleapoutput:
            line = line.decode('utf-8')
            r = re.search("^R<(.*) (.*)>", line)
            if r:
                if r.group(1) not in reslist:
                    reslist[r.group(1)] = []
                reslist[r.group(1)].append(r.group(2))
        return reslist
    elif returnlist == 'WAT':
        return reslist
    else:
        return wateradded


def solvate(tleapfile, pdbfile=None, pbctype=1, bufferwater='12.0A',
            waterbox='TIP3PBOX', neutralize=1, counter_cation='Na+', counter_anion='Cl-',
            addions=None, saveprefix='solvated'):
    """
    This routine solvates a solute system with a specified amount of water/ions.

    ARGUMENTS

    tleapfile : a fully functioning tleap file which is capable of preparing the system in gas phase. It should load all parameter files that will be necessary for solvation, including the water model and any custom residues. Assumes the final conformations of the solute are loaded via PDB.

    pdbfile : if present, the function will search for any loadpdb commands in the tleapfile and replace whatever is there with pdbfile.  This would be used for cases where we have modified PDBs.

    pbctype : the type of periodic boundary conditions. 0 = cubic, 1 = rectangular, 2 = truncated octahedron.  If rectangular, only the z-axis buffer can be manipulated; the x- and y-axis will use a 10 Ang buffer.

    waterbox : the water box name to use with the solvatebox/solvateoct command. 

    neutralize : 0 = do not neutralize the system, 1 = neutralize the system. the counterions to be used are specified below with 'counter_catio' and 'counter_anion'.

    counter_cation : a mask to specify neutralizing cations

    counter_anion : a mask to specify neturalizing anions

    addions : a list of residue masks and values which indicate how much additional ions to add. The format for the values is as following: if the value is an integer, then add that exact integer amount of ions; if the value is followed by an 'M', then add that amount in molarity;  if 'm', add by molality.  example: ['Mg+',5, 'Cl-',10, 'K+','0.050M']

    """

    unit = 'model'
    addion_residues = []
    addion_values = []
    bufferval = [0.0]
    buffer_iterexp = 1
    wat_added = [0.0]

    dir = os.path.dirname(tleapfile)
    if dir == '':
        dir = '.'

    # Read tleapfile
    tleaplines = []
    with open(tleapfile, 'r') as f:
        for line in f.readlines():
            if re.search('loadpdb', line):
                words = line.rstrip().replace('=', ' ').split()
                if pdbfile is None:
                    pdbfile = words[2]
                    log.debug('PDB: {}'.format(pdbfile))
                unit = words[0]
                tleaplines.append(
                    "{} = loadpdb {}\n".format(unit, pdbfile))
            if not re.search(r"^\s*addions|^\s*addions2|^\s*addionsrand|^\s*desc|"
                             r"^\s*quit|^\s*solvate|loadpdb|^\s*save", line, re.IGNORECASE):
                tleaplines.append(line)

    # If bufferwater ends with 'A', meaning it is a buffer distance...
    if str(bufferwater).endswith('A'):
        # Let's get a rough value of the number of waters if the buffer value
        # is given as a string.
        bufferval.append(float(bufferwater[:-1]))
        write_tleapin(filename='tleap_apr_solvate.in',
                      directory=dir,
                      tleaplines=tleaplines,
                      unit=unit,
                      pbctype=pbctype,
                      bufferval=bufferval[-1],
                      waterbox=waterbox,
                      neutralize=0,
                      counter_cation=None,
                      counter_anion=None,
                      addion_residues=None,
                      addion_values=None,
                      removewat=None,
                      saveprefix=None
                     )
        log.debug('Getting an itial estimate for how many waters fit in {}...'.format(bufferwater))
        bufferwater = countresidues(filename='tleap_apr_solvate.in', directory=dir)

    if addions:
        if len(addions) % 2 == 1:
            raise Exception("Error: The 'addions' list requires an even number of elements. "
                            "Make sure there is a residue mask followed by a value for "
                            "each ion to be added")
        for i, txt in enumerate(addions):
            if i % 2 == 0:
                addion_residues.append(txt)
            else:
                # Temporary treat m and M identical
                if str(txt).endswith('m') or str(txt).endswith('M'):
                    # Figure out number of ions for desired molality
                    addion_values.append(
                        float(txt[:-1]) * float(bufferwater) * 0.018)
                else:
                    addion_values.append(txt)

    # First adjust wat_added by changing the bufferval
    cycle = 0
    buffer_iter = 0
    while cycle < 50:
        write_tleapin(filename='tleap_apr_solvate.in',
                      directory=dir,
                      tleaplines=tleaplines,
                      unit=unit,
                      pbctype=pbctype,
                      bufferval=bufferval[-1],
                      waterbox=waterbox,
                      neutralize=neutralize,
                      counter_cation=counter_cation,
                      counter_anion=counter_anion,
                      addion_residues=addion_residues,
                      addion_values=addion_values,
                      removewat=None,
                      saveprefix=None
                      )
        wat_added.append(countresidues(filename='tleap_apr_solvate.in', directory=dir))
        print(cycle,buffer_iter,":",bufferwater,':',bufferval[-1],buffer_iterexp,wat_added[-2],wat_added[-1])
        cycle += 1
        buffer_iter += 1
        if 0 <= (wat_added[-1] - bufferwater) < 12 or buffer_iterexp < -3:
            # Possible failure mode: if the tolerance here is very small (0 < () < 1),
            # the loop can exit with bufferval that adds fewer waters than
            # bufferwater
            log.info('Done!')
            break
        # Possible location of a switch to adjust the bufferval by polynomial
        # fit approach.
        elif wat_added[-2] > bufferwater and wat_added[-1] > bufferwater:
            bufferval.append(bufferval[-1] + -1 * (10**buffer_iterexp))
        elif wat_added[-2] > bufferwater and wat_added[-1] < bufferwater:
            if buffer_iter > 1:
                buffer_iterexp -= 1
                buffer_iter = 0
                bufferval.append(bufferval[-1] + 5 * (10**buffer_iterexp))
            else:
                bufferval.append(bufferval[-1] + 1 * (10**buffer_iterexp))
        elif wat_added[-2] < bufferwater and wat_added[-1] > bufferwater:
            if buffer_iter > 1:
                buffer_iterexp -= 1
                buffer_iter = 0
                bufferval.append(bufferval[-1] + -5 * (10**buffer_iterexp))
            else:
                bufferval.append(bufferval[-1] + -1 * (10**buffer_iterexp))
        elif wat_added[-2] < bufferwater and wat_added[-1] < bufferwater:
            bufferval.append(bufferval[-1] + 1 * (10**buffer_iterexp))
        else:
            raise Exception(
                "The bufferval search died due to an unanticipated set of variable values")

    if cycle >= 50:
        raise Exception(
            "Automatic adjustment of the buffer value was unable to converge on \
            a solution with sufficient tolerance")
    elif wat_added[-1] - bufferwater < 0:
        raise Exception(
            "Automatic adjustment of the buffer value resulted in fewer waters \
            added than targeted by bufferwater. Try increasing the tolerance in the above loop")
    else:
        watover = 0
        cycle = 0
        while wat_added[-1] != bufferwater or cycle == 0:
            # Note I don't think there should be water removal errors, but if
            # so, this loop and '+=' method is an attempt to fix.
            watover += wat_added[-1] - bufferwater
            watlist = countresidues(
                filename='tleap_apr_solvate.in', returnlist='WAT', directory=dir)
            if watover == 0:
                removewat = None
            else:
                removewat = watlist[-1 * watover:]

            write_tleapin(filename='tleap_apr_solvate.in',
                          directory=dir,
                          tleaplines=tleaplines,
                          unit=unit,
                          pbctype=pbctype,
                          bufferval=bufferval[-1],
                          waterbox=waterbox,
                          neutralize=neutralize,
                          counter_cation=counter_cation,
                          counter_anion=counter_anion,
                          addion_residues=addion_residues,
                          addion_values=addion_values,
                          removewat=removewat,
                          saveprefix=saveprefix
                         )
            reslist = countresidues(filename='tleap_apr_solvate.in',
                                    returnlist='ALL', directory=dir)
            wat_added.append(len(reslist['WAT']))
            for key, value in sorted(reslist.items()):
                log.info('{}\t{}'.format(key, len(value)))
            cycle += 1
            if cycle >= 10:
                raise Exception(
                    "Solvation failed due to an unanticipated problem with water removal")


# solvate('../test/cb6-but/tleap.in', bufferwater=2003, pbctype=1, addions=['K+', 5, 'BR', 5])
