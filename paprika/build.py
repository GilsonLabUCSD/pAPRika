import os as os
import re as re
import subprocess as sp

import logging as log
import numpy as np
import parmed as pmd
from parmed.structure import Structure as ParmedStructureClass
from paprika import utils

N_A = 6.0221409 * 10**23
angstrom_cubed_to_liters = 1 * 10**-27


def read_tleaplines(pdb_file=None, path='./', filename='tleap.in', filepath=None):
    """
    Read a tleap input file and return a list containing each line of instruction.
    
    Parameters:
    ----------
    pdb_file : {str}, optional
        The file name of a to-be-processed `pdb` file, otherwise detected from the input file.
    path : {str}, optional
        The directory of the output file, if `filepath` is not specified (the default is './')
    filename : {str}, optional
        The name of the output file, if `filepath` is not specified (the default is 'dummy.mol2')
    filepath : {str}, optional
        The full path (directory and file) of the output (the default is None, which means `path` and `filename` will be used)

    Returns:
    -------
    lines : {list}
        The list of lines in a `tleap` input file
    """

    if filepath is None:
        filepath = path + filename

    lines = []
    with open(filepath, 'r') as f:
        for line in f.readlines():
            if re.search('loadpdb', line):
                words = line.rstrip().replace('=', ' ').split()
                if pdb_file is None:
                    pdb_file = words[2]
                unit = words[0]
                lines.append("{} = loadpdb {}\n".format(unit, pdb_file))
            if not re.search(r"^\s*addions|^\s*addions2|^\s*addionsrand|^\s*desc|"
                             r"^\s*quit|^\s*solvate|loadpdb|^\s*save", line, re.IGNORECASE):
                lines.append(line)

    return lines


def write_tleapin(lines, options, path='./', filename='tleap.in', filepath=None):
    """

    Parameters
    ----------
    lines : list
        Boilerplate `tleap` input file passed as a list of lines.
    unit : str
        The "model" during processing (see the AMBER manual for more information)
    pbc_type : int
        The solvation model passed along to `tleap`
    buffer_values : float
        The amount of solvation buffer
    water_box : str
        The solvation model used (see the AMBER manual for allowed values)
    neutralize : bool
        Whether to enforce charge neutrality in the solvated structure
    counter_cation : str
        The positive ion used for charge neutrality (see the AMBER manual for allowed values)
    counter_anion : str
        The negative ion used for charge neutrality (see the AMBER manual for allowed values)
    addion_residues : str
        Any additional ions added to the solvated structure
    addion_values : list
        The amount of any additional ions added to the solvated structure
    remove_wat : None or list
        If none: do nothing. If list: a list of water residues to remove from the solvated structure.
    output_prefix : str
        The file name of the solvated `prmtop` and `rst7` files
    path : {str}, optional
        The directory of the output file, if `filepath` is not specified (the default is './')
    filename : {str}, optional
        The name of the output file, if `filepath` is not specified (the default is 'tleap.in')
    filepath : {str}, optional
        The full path (directory and file) of the output (the default is None, which means `path` and `filename` will be used)

    """

    # Kludge to get things working now...
    unit = options['unit']
    pbc_type = options['pbc_type']
    buffer_values = options['buffer_values']
    water_box = options['water_box']
    neutralize = options['neutralize']
    counter_cation = options['counter_cation']
    counter_anion = options['counter_anion']
    addion_residues = options['addion_residues']
    addion_values = options['addion_values']
    remove_water = options['remove_water']
    output_prefix = options['output_prefix']

    if filepath is None:
        filepath = path + filename

    with open(filepath, 'w') as f:
        for line in lines:
            f.write(line)
        if pbc_type == 0:
            f.write("solvatebox {} {} {:0.5f} iso\n".format(unit, water_box, buffer_values))
        elif pbc_type == 1:
            f.write("solvatebox {} {} {{10.0 10.0 {:0.5f}}}\n".format(unit, water_box, buffer_values))
        elif pbc_type == 2:
            f.write("solvateoct {} {} {:0.5f} iso\n".format(unit, water_box, buffer_values))
        elif pbc_type is None:
            f.write("# Skipping solvation ...\n")
        else:
            raise Exception(
                "Incorrect pbctype value provided: " + str(pbc_type) + ". Only 0, 1, 2, and None are valid")
        if neutralize:
            f.write("addionsrand {} {} 0\n".format(unit, counter_cation))
            f.write("addionsrand {} {} 0\n".format(unit, counter_anion))
        if isinstance(addion_residues, list):
            try:
                assert len(addion_residues) == len(addion_values)
            except AssertionError as e:
                raise Exception("Different number of additional ion residues and amounts to add to the structure.")
            for i, res in enumerate(addion_residues):
                f.write("addionsrand {} {} {}\n".format(unit, res, addion_values[i]))
        if remove_wat:
            for water_number in remove_wat:
                f.write("remove {} {}.{}\n".format(unit, unit, water_number))

        f.write("savepdb {} {}.pdb\n".format(unit, output_prefix))
        f.write("saveamberparm {} {}.prmtop {}.rst7\n".format(unit, output_prefix, output_prefix))
        f.write("desc {}\n".format(unit))
        f.write("quit\n")


def run_tleap(path='./', filename='tleap.in'):
    """
    Execute tLEaP.
    """

    utils.check_for_leap_log(path=path)
    p = sp.Popen('tleap -s -f ' + filename, stdout=sp.PIPE, bufsize=1, universal_newlines=True, cwd=path, shell=True)
    output = []
    while p.poll() is None:
        line = p.communicate()[0]
        output.append(line)
    if p.poll() is None:
        p.kill()

    return output


def basic_tleap(input_file='tleap.in', input_path='./', output_prefix='solvate', output_path='./', pdb_file=None):
    """
    Run tleap with a user supplied tleap script, optionally substitute PDB.
    """

    log.debug('Reading {}/{}, writing {}/{}.in, and executing ...'.format(input_path, input_file, output_path,
                                                                          output_prefix))

    lines = read_tleaplines(pdb_file=pdb_file, path=input_path, filename=input_file)
    options = {}
    options['pbc_type'] = None
    options['neutralize'] = False
    options['output_prefix'] = output_prefix
    write_tleapin(lines, options, path=output_path, filename=output_prefix + '.in')
    run_tleap(path=output_path, filename=output_prefix + '.in')


def count_residues(path='./', filename='tleap.in'):
    """Run and parse `tleap` output and return a dictionary of added residues.
    
    """
    output = run_tleap(path=path, filename=filename)
    # Reurn a dictionary of {'RES' : number of RES}
    residues = {}
    for line in output[0].splitlines():
        # Is this line a residue from `desc` command?
        match = re.search("^R<(.*) ", line)
        if match:
            residue_name = match.group(1)
            # If this residue is not in the dictionary, initialize and
            # set the count to 1.
            if residue_name not in residues:
                residues[residue_name] = 1
            # If this residue is in the dictionary, increment the count
            # each time we find an instance.
            elif residue_name in residues:
                residues[residue_name] += 1
    return residues


def count_volume(path='./', filename='tleap.in'):
    output = run_tleap(path=path, filename=filename)
    # Return the total simulation volume
    for line in output[0].splitlines():
        match = re.search("Volume(.*)>", line)
        if match:
            volume = float(r.group(1)[1:-4])
            return volume
    log.warning('Could not determine total simulation volume.')
    return None


def count_waters(path='./', filename='tleap.in'):
    """
    Run and parse `tleap` output and return the number of water residues.
    """
    output = run_tleap(path=path, filename=filename)

    # Return a list of residue numbers for the waters
    water_residues = []
    for line in output[0].splitlines():
        # Is this line a water?
        match = re.search("^R<WAT (.*)>", line)
        if match:
            water_residues.append(match.group(1))
    return water_residues


def solvate(tleap_file,
            pdb_file=None,
            pbctype=1,
            bufferwater='12.0A',
            waterbox='TIP3PBOX',
            neutralize=True,
            counter_cation='Na+',
            counter_anion='Cl-',
            addions=None,
            saveprefix='solvated',
            path='./'):
    """
    This routine solvates a solute system with a specified amount of water/ions.

    ARGUMENTS

    tleapfile : a fully functioning tleap file which is capable of preparing
    the system in gas phase. It should load all parameter files that will be
    necessary for solvation, including the water model and any custom
    residues. Assumes the final conformations of the solute are loaded via
    PDB.

    pdb_file : if present, the function will search for any loadpdb commands in
    the tleapfile and replace whatever is there with pdb_file.  This would be
    used for cases where we have modified PDBs. returnlist can be ALL for all
    residuse or WAT for waters.

    pbctype : the type of periodic boundary conditions. 0 = cubic, 1 =
    rectangular, 2 = truncated octahedron, None = no solvation.  If
    rectangular, only the z-axis buffer can be manipulated; the x- and y-axis
    will use a 10 Ang buffer.

    waterbox : the water box name to use with the solvatebox/solvateoct
    command.

    neutralize : False = do not neutralize the system, True = neutralize the
    system. the counterions to be used are specified below with
    'counter_cation' and 'counter_anion'.

    counter_cation : a mask to specify neutralizing cations

    counter_anion : a mask to specify neturalizing anions

    addions : a list of residue masks and values which indicate how much
    additional ions to add. The format for the values is as following: if the
    value is an integer, then add that exact integer amount of ions; if the
    value is followed by an 'M', then add that amount in molarity;  if 'm',
    add by molality.  example: ['Mg+',5, 'Cl-',10, 'K+','0.050M']

    """

    unit = 'model'
    addion_residues = []
    addion_values = []
    buffer_values = [0.0]
    buffer_iterexp = 1
    wat_added = [0.0]

    # Read tleapfile
    lines = read_tleaplines(pdb_file=pdb_file, path=path, filename=tleap_file)

    # If bufferwater ends with 'A', meaning it is a buffer distance...
    # NMH: there is no other branch to this `if` statement...
    if str(bufferwater).endswith('A'):
        # Let's get a rough value of the number of waters if the buffer value
        # is given as a string.
        buffer_values.append(float(bufferwater[:-1]))

        options = {}
        options['pbc_type'] = pbc_type
        options['buffer_values'] = buffer_values[-1]
        options['water_box'] = waterbox

        write_tleapin(
            filename='tleap_apr_solvate.in',
            path=path,
            tleaplines=tleaplines,
            unit=unit,
            pbc_type=pbctype,
            buffer_values=buffer_values[-1],
            water_box=waterbox,
            neutralize=False,
            counter_cation=None,
            counter_anion=None,
            addion_residues=None,
            addion_values=None,
            remove_wat=None,
            output_prefix=None)
        # Get the number of water residues added...
        bufferwater = count_residues(
            filename='tleap_apr_solvate.in', directory=path, choice='residue_dictionary')['WAT']
    if addions:
        if len(addions) % 2 == 1:
            raise Exception("Error: The 'addions' list requires an even number of elements. "
                            "Make sure there is a residue mask followed by a value for "
                            "each ion to be added")
        for i, txt in enumerate(addions):
            if i % 2 == 0:
                addion_residues.append(txt)
            else:
                # User specifies molaliy...
                if str(txt).endswith('m'):
                    # number to add = (molality) x (number waters) x (0.018 kg/mol per water)
                    addion_values.append(int(np.ceil(float(txt[:-1]) * float(bufferwater) * 0.018)))
                # User specifies molarity...
                elif str(txt).endswith('M'):
                    volume = count_residues(
                        filename='tleap_apr_solvate.in', directory=path, choice='residue_dictionary', volume=True)
                    number_of_atoms = float(txt[:-1]) * N_A
                    liters = volume * angstrom_cubed_to_liters
                    atoms_to_add = number_of_atoms * liters
                    addion_values.append(np.ceil(atoms_to_add))
                else:
                    addion_values.append(int(txt))

    # First adjust wat_added by changing the buffer_values
    cycle = 0
    buffer_iter = 0
    while cycle < 50:
        write_tleapin(
            filename='tleap_apr_solvate.in',
            directory=path,
            tleaplines=tleaplines,
            unit=unit,
            pbc_type=pbctype,
            buffer_values=buffer_values[-1],
            water_box=waterbox,
            neutralize=neutralize,
            counter_cation=counter_cation,
            counter_anion=counter_anion,
            addion_residues=addion_residues,
            addion_values=addion_values,
            remove_wat=None,
            output_prefix=None)
        # Get the number of water residues added...
        wat_added.append(
            count_residues(filename='tleap_apr_solvate.in', directory=path, choice='residue_dictionary')['WAT'])
        log.debug(cycle, buffer_iter, ":", bufferwater, ':', buffer_values[-1], buffer_iterexp, wat_added[-2],
                  wat_added[-1])
        cycle += 1
        buffer_iter += 1
        if 0 <= (wat_added[-1] - bufferwater) < 12 or \
                (buffer_iterexp < -3 and (wat_added[-1] - bufferwater) > 0):
            # Possible failure mode: if the tolerance here is very small (0 < () < 1),
            # the loop can exit with buffer_values that adds fewer waters than
            # bufferwater
            log.info('Done!')
            break
        # Possible location of a switch to adjust the buffer_values by polynomial
        # fit approach.
        elif wat_added[-2] > bufferwater and wat_added[-1] > bufferwater:
            buffer_values.append(buffer_values[-1] + -1 * (10**buffer_iterexp))
        elif wat_added[-2] > bufferwater and wat_added[-1] < bufferwater:
            if buffer_iter > 1:
                buffer_iterexp -= 1
                buffer_iter = 0
                buffer_values.append(buffer_values[-1] + 5 * (10**buffer_iterexp))
            else:
                buffer_values.append(buffer_values[-1] + 1 * (10**buffer_iterexp))
        elif wat_added[-2] < bufferwater and wat_added[-1] > bufferwater:
            if buffer_iter > 1:
                buffer_iterexp -= 1
                buffer_iter = 0
                buffer_values.append(buffer_values[-1] + -5 * (10**buffer_iterexp))
            else:
                buffer_values.append(buffer_values[-1] + -1 * (10**buffer_iterexp))
        elif wat_added[-2] < bufferwater and wat_added[-1] < bufferwater:
            buffer_values.append(buffer_values[-1] + 1 * (10**buffer_iterexp))
        else:
            raise Exception("The buffer_values search died due to an unanticipated set of variable values")

    if cycle >= 50:
        raise Exception("Automatic adjustment of the buffer value was unable to converge on \
            a solution with sufficient tolerance")
    elif wat_added[-1] - bufferwater < 0:
        raise Exception("Automatic adjustment of the buffer value resulted in fewer waters \
            added than targeted by bufferwater. Try increasing the tolerance in the above loop")
    else:
        watover = 0
        cycle = 0
        while wat_added[-1] != bufferwater or cycle == 0:
            # Note I don't think there should be water removal errors, but if
            # so, this loop and '+=' method is an attempt to fix.
            watover += wat_added[-1] - bufferwater
            watlist = count_residues(filename='tleap_apr_solvate.in', choice='water_residues', directory=path)
            if watover == 0:
                removewat = None
            else:
                removewat = watlist[-1 * watover:]

            write_tleapin(
                filename='tleap_apr_solvate.in',
                directory=path,
                tleaplines=tleaplines,
                unit=unit,
                pbc_type=pbctype,
                buffer_values=buffer_values[-1],
                water_box=waterbox,
                neutralize=neutralize,
                counter_cation=counter_cation,
                counter_anion=counter_anion,
                addion_residues=addion_residues,
                addion_values=addion_values,
                remove_wat=removewat,
                output_prefix=saveprefix)
            # Get the full residue dictionary...
            reslist = count_residues(filename='tleap_apr_solvate.in', choice='residue_dictionary', directory=path)
            wat_added.append(reslist['WAT'])
            for key, value in sorted(reslist.items()):
                log.info('{}\t{}'.format(key, value))
            cycle += 1
            if cycle >= 10:
                raise Exception("Solvation failed due to an unanticipated problem with water removal")
