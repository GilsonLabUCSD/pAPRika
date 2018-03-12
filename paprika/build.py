import os as os
import re as re
import subprocess as sp

import logging as log
import numpy as np
import parmed as pmd
from parmed.structure import Structure as ParmedStructureClass
from paprika import utils

N_A = 6.0221409 * 10**23
ANGSTROM_CUBED_TO_LITERS = 1 * 10**-27


def default_tleap_options():
    """
    Return a dictionary of default `tleap` options.
    Returns
    -------

    """
    options = {}
    options['unit'] = 'model'
    options['pbc_type'] = 1
    options['buffer_value'] = []
    options['water_box'] = 'TIP3PBOX'
    options['neutralize'] = False
    options['counter_cation'] = None
    options['counter_anion'] = None
    options['add_ion_residues'] = None
    options['remove_water'] = None
    options['output_prefix'] = 'solvate'
    options['path'] = './'
    return options


def read_tleap_lines(pdb_file=None, path='./', file_name='tleap.in', file_path=None):
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

    if file_path is None:
        file_path = path + file_name

    lines = []
    with open(file_path, 'r') as f:
        for line in f.readlines():
            if re.search('loadpdb', line):
                words = line.rstrip().replace('=', ' ').split()
                if pdb_file is None:
                    pdb_file = words[2]
                unit = words[0]
                lines.append("{} = loadpdb {}\n".format(unit, pdb_file))
            # Skip over any included solvation and ionization commands...
            if not re.search(r"^\s*addions|^\s*addions2|^\s*addionsrand|^\s*desc|"
                             r"^\s*quit|^\s*solvate|loadpdb|^\s*save", line, re.IGNORECASE):
                lines.append(line)

    return lines


def write_tleapin(lines, options):
    """

    Parameters
    ----------
    lines : list
        Boilerplate `tleap` input file passed as a list of lines.
    options : dict
        Dictionary containing `tleap` options, set up by `default_tleap_options()`
    """

    file_path = options['path'] + options['output_prefix'] + '.in'
    with open(file_path, 'w') as f:
        for line in lines:
            f.write(line)

        if options['pbc_type'] == 0:
            f.write("solvatebox {} {} {} iso\n".format(options['unit'], options['water_box'],
                                                            options['buffer_value']))
        elif options['pbc_type'] == 1:
            f.write("solvatebox {} {} {{10.0 10.0 {}}}\n".format(options['unit'], options['water_box'],
                                                                      options['buffer_value']))
        elif options['pbc_type'] == 2:
            f.write("solvateoct {} {} {} iso\n".format(options['unit'], options['water_box'],
                                                            options['buffer_value']))
        elif options['pbc_type'] is None:
            f.write("# Skipping solvation ...\n")
        else:
            raise Exception(
                "Incorrect pbctype value provided: " + str(options['pbc_type']) + ". Only 0, 1, 2, and None are valid")
        if options['neutralize']:
            f.write("addionsrand {} {} 0\n".format(options['unit'], options['counter_cation']))
            f.write("addionsrand {} {} 0\n".format(options['unit'], options['counter_anion']))
        # Additional ions should be specified as a list, with residue name and number of ions in pairs, like ['Na+',
        # 5] for five additional sodium ions. By this point, if the user specified a molality or molarity,
        # it should already have been converted into a number.
        if options['add_ion_residues']:
            for residue, amount in zip(options['add_ion_residues'][0::2],
                                       options['add_ion_residues'][1::2]):
                f.write("addionsrand {} {} {}\n".format(options['unit'], residue, amount))
        if options['remove_water']:
            for water_number in options['remove_water']:
                f.write("remove {} {}.{}\n".format(options['unit'], options['unit'], water_number))

        if options['output_prefix']:
            f.write("savepdb {} {}.pdb\n".format(options['unit'], options['output_prefix']))
            f.write("saveamberparm {} {}.prmtop {}.rst7\n".format(options['unit'], options['output_prefix'],
                                                                  options['output_prefix']))
        else:
            pass
        f.write("desc {}\n".format(options['unit']))
        f.write("quit\n")


def run_tleap(path='./', file_name='tleap.in'):
    """
    Execute tLEaP.
    """

    utils.check_for_leap_log(path=path)
    p = sp.Popen('tleap -s -f ' + file_name, stdout=sp.PIPE, bufsize=1, universal_newlines=True, cwd=path, shell=True)
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

    lines = read_tleap_lines(pdb_file=pdb_file, path=input_path, file_name=input_file)
    options = default_tleap_options()
    options['path'] = output_path
    options['output_prefix'] = output_prefix
    write_tleapin(lines, options)
    run_tleap(path=output_path, file_name=output_prefix + '.in')


def count_residues(path='./', file_name='tleap.in'):
    """Run and parse `tleap` output and return a dictionary of added residues.
    
    """
    output = run_tleap(path=path, file_name=file_name)
    # Return a dictionary of {'RES' : number of RES}
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


def count_volume(path='./', file_name='tleap.in'):
    output = run_tleap(path=path, file_name=file_name)
    # Return the total simulation volume
    for line in output[0].splitlines():
        line = line.strip()
        if "Volume" in line:
            match = re.search("Volume(.*)", line)
            volume = float(match.group(1)[1:-4])
            return volume
    log.warning('Could not determine total simulation volume.')
    return None


def count_waters(path='./', file_name='tleap.in'):
    """
    Run and parse `tleap` output and return the number of water residues.
    """
    output = run_tleap(path=path, file_name=file_name)

    # Return a list of residue numbers for the waters
    water_residues = []
    for line in output[0].splitlines():
        # Is this line a water?
        match = re.search("^R<WAT (.*)>", line)
        if match:
            water_residues.append(match.group(1))
    return water_residues


def quick_check(lines, options):
    """
    Quickly check the number of waters added for a given buffer size.
    """
    write_tleapin(lines, options)
    waters = count_residues(file_name=options['output_prefix'] + '.in', path=options['path'])['WAT']
    return waters

def set_target_number_of_waters(lines, options, buffer_target):
    # If buffer_water ends with 'A', meaning it is a buffer distance...
    if isinstance(buffer_target, str) and buffer_target[-1] == 'A':
        # Let's get a rough value of the number of waters if the buffer target is given as a string.
        options['buffer_value'] = buffer_target[:-1]
        waters = quick_check(lines, options)
        log.debug('Initial guess of {} waters for a buffer size of {}...'.format(waters, buffer_target))
        # This is now the target number of waters for solvation...
        return waters
    elif isinstance(buffer_target, int):
        return buffer_target
    # Otherwise, the number of waters to add is specified as an integer, not a distance...
    else:
        raise Exception("The `buffer_target` should either be a string ending with 'A' (e.g., 12A) for 12 Angstroms of "
                        "buffer or an int (e.g., 2000) for 2000 waters.")


def set_additional_ions(add_ions, options, buffer_target):
    if not add_ions:
        return None
    if len(add_ions) < 2:
        raise Exception("No amount specified for additional ions.")
    if len(add_ions) % 2 == 1:
        raise Exception("The 'add_ions' list requires an even number of elements. "
                        "Make sure there is a residue mask followed by a value for "
                        "each ion to be added (or molarity ending in 'M' or molality ending in 'm').")
    add_ion_residues = []
    for ion, amount in zip(add_ions[0::2], add_ions[1::2]):
        add_ion_residues.append(ion)
        if isinstance(amount, int):
            add_ion_residues.append(amount)
        elif isinstance(amount, str) and amount[-1] == 'm':
            # User specifies molality...
            # number to add = (molality) x (number waters) x (0.018 kg/mol per water)
            number_to_add = int(np.ceil(float(amount[:-1]) * buffer_target * 0.018))
            add_ion_residues.append(number_to_add)
        elif isinstance(amount, str) and amount[-1] == 'M':
            # User specifies molarity...
            volume = count_volume(file_name=options['output_prefix'] + '.in', path=options['path'])
            number_of_atoms = float(amount[:-1]) * N_A
            liters = volume * ANGSTROM_CUBED_TO_LITERS
            number_to_add = number_of_atoms * liters
            add_ion_residues.append(number_to_add)
    return add_ion_residues


def adjust_buffer_value(number_of_waters, target_number_of_waters, buffer_values, exponent):

        # If the last two rounds of solvation have too many waters, make the buffer smaller...
        if number_of_waters[-2] > target_number_of_waters and number_of_waters[-1] > target_number_of_waters:
            log.debug('Adjustment loop 1')
            return buffer_values[-1] + -1 * (10**exponent), exponent

        # If the last two rounds of solvation had too few waters, make the buffer bigger...
        elif number_of_waters[-2] < target_number_of_waters and number_of_waters[-1] < target_number_of_waters:
            log.debug('Adjustment loop 2')
            return buffer_values[-1] + 1 * (10**exponent), exponent

        # If the number of waters was greater than the target and is now less than the target, make the buffer a bit
        # bigger, by an increasingly smaller amount...
        elif number_of_waters[-2] > target_number_of_waters and number_of_waters[-1] < target_number_of_waters:
            log.debug('Adjustment loop 3')
            exponent -= 1
            return buffer_values[-1] + 5 * (10**exponent), exponent

        # If the number of waters was less than the target and is now greater than the target, make the buffer a bit
        # smaller, by an increasingly bigger amount...
        elif number_of_waters[-2] < target_number_of_waters and number_of_waters[-1] > target_number_of_waters:
            log.debug('Adjustment loop 4')
            exponent -= 1
            return buffer_values[-1] + -5 * (10**exponent), exponent
        else:
            raise Exception("The buffer_values search died due to an unanticipated set of variable values")


def remove_waters_manually(lines, number_of_waters, target_number_of_waters, options):
    cycle = 0
    max_cycles = 10
    waters = number_of_waters[-1]
    while waters > target_number_of_waters:
        water_surplus = (waters - target_number_of_waters)
        water_residues = count_waters(file_name=options['output_prefix'] + '.in', path=options['path'])
        waters_to_remove = water_residues[-1 * water_surplus:]
        log.debug('Manually removing waters... {}'.format(waters_to_remove))
        options['remove_water'] = waters_to_remove
        write_tleapin(lines, options)

        residues = count_residues(file_name=options['output_prefix'] + '.in', path=options['path'])
        waters = residues['WAT']
        if waters == target_number_of_waters:
            for key, value in sorted(residues.items()):
                log.info('{}\t{}'.format(key, value))
            return
        cycle += 1
        if cycle > max_cycles:
            raise Exception("Solvation failed due to an unanticipated problem with water removal.")


def solvate(tleap_file,
            pdb_file=None,
            pbc_type=1,
            buffer_target='12.0A',
            water_box='TIP3PBOX',
            neutralize=True,
            counter_cation='Na+',
            counter_anion='Cl-',
            add_ions=None,
            output_prefix='solvate',
            path='./'):
    """

    """

    # Read template and setup default `tleap` options.
    lines = read_tleap_lines(pdb_file=pdb_file, path=path, file_name=tleap_file)
    options = default_tleap_options()
    options['pbc_type'] = pbc_type
    options['pdb_file'] = pdb_file
    options['neutralize'] = neutralize
    options['counter_cation'] = counter_cation
    options['counter_anion'] = counter_anion
    options['output_prefix'] = output_prefix
    options['water_box'] = water_box
    options['path'] = path
    # If `buffer_target` is a string ending with 'A', an estimate of the number of waters is generated, otherwise,
    # the target is returned.
    target_number_of_waters = set_target_number_of_waters(lines, options, buffer_target)

    if add_ions:
        options['add_ion_residues'] = set_additional_ions(add_ions, options, target_number_of_waters)

    # First, a coarse adjustment...
    # This will run for 50 iterations or until we (a) have more waters than the target and (b) are within ~12 waters
    # of the target (that can be manually removed).
    cycle = 0
    max_cycles = 50
    # Start small (I think this will work)...
    options['buffer_value'] = 1
    # List of buffer values attempted...
    buffer_values = [0]
    # Initial exponent used to narrow the search of buffer values to give a target number of waters...
    exponent = 1
    # Number of waters corresponding to each buffer value...
    number_of_waters = [0]
    while cycle < max_cycles:

        # Start by not manually removing any water, just adjusting the buffer value to get near the target number of
        # waters...
        options['remove_water'] = None
        write_tleapin(lines, options)
        # Find out how many waters for *this* buffer value...
        residues = count_residues(file_name=output_prefix + '.in', path=path)
        log.debug(residues)
        waters = residues['WAT']
        number_of_waters.append(waters)
        buffer_values.append(options['buffer_value'])

        log.debug('Cycle {0:02d}\t'.format(cycle), options['buffer_value'], waters,
                  '({})'.format(target_number_of_waters))

        # Possible location of a switch to adjust the buffer_values by polynomial
        # fit approach.

        # If we've nailed it, break!
        if waters == target_number_of_waters:
            return
        # If we are close, go to fine adjustment...
        elif waters > target_number_of_waters and (waters - target_number_of_waters) < 12:
            remove_waters_manually(lines, number_of_waters, target_number_of_waters, options)
            return
        # Otherwise, try to keep adjusting the number of waters...
        else:
            options['buffer_value'], exponent = adjust_buffer_value(number_of_waters, target_number_of_waters, buffer_values,
                                                         exponent)
            cycle += 1


    if cycle >= max_cycles and waters > target_number_of_waters:
        remove_waters_manually(lines, number_of_waters, target_number_of_waters, options)

    if cycle >= max_cycles and waters < target_number_of_waters:
        raise Exception("Automatic adjustment of the buffer value resulted in fewer waters \
            added than targeted by `buffer_water`. Try increasing the tolerance in the above loop")
    else:
        raise Exception("Automatic adjustment of the buffer value was unable to converge on \
            a solution with sufficient tolerance")