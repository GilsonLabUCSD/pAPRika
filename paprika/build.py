import os as os
import re as re
import subprocess as sp

import logging as log
import numpy as np
import parmed as pmd
from parmed.structure import Structure as ParmedStructureClass
from paprika import utils

N_A = 6.0221409 * 10 ** 23
ANGSTROM_CUBED_TO_LITERS = 1 * 10 ** -27


def default_tleap_options():
    """
    Return a dictionary of default `tleap` options.

    """
    options = {}
    options['unit'] = 'model'
    options['pbc_type'] = 'rectangular'
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
    Read a `tleap` input file and return a list containing each line of instruction, minus solvation, which we re-write
    later.
    
    Parameters:
    ----------
    pdb_file : {str}, optional
        The file name of a to-be-processed `pdb` file, otherwise detected from the input file.
    path : {str}, optional
        The directory of the output file, if `filepath` is not specified (the default is './')
    file_name : {str}, optional
        The name of the output file, if `filepath` is not specified (the default is 'dummy.mol2')
    file_path : {str}, optional
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
    Write a `tleap` input file using lines from a template file and a dictionary of options.

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

        if 'cubic' in options['pbc_type']:
            f.write("solvatebox {} {} {} iso\n".format(options['unit'], options['water_box'],
                                                       options['buffer_value']))
        elif 'rectangular' in options['pbc_type']:
            f.write("solvatebox {} {} {{10.0 10.0 {}}}\n".format(options['unit'], options['water_box'],
                                                                 options['buffer_value']))
        elif 'octahedral' in options['pbc_type']:
            f.write("solvateoct {} {} {} iso\n".format(options['unit'], options['water_box'],
                                                       options['buffer_value']))
        elif options['pbc_type'] is None:
            f.write("# Skipping solvation ...\n")
        else:
            raise Exception(
                "Incorrect pbctype value provided: " + str(options['pbc_type']) + ". Only `cubic`, `rectangular`, "
                                                                                  "`octahedral`, and "
                                                                                  "None are valid")
        if options['neutralize']:
            f.write("addionsrand {} {} 0\n".format(options['unit'], options['counter_cation']))
            f.write("addionsrand {} {} 0\n".format(options['unit'], options['counter_anion']))
        # Additional ions should be specified as a list, with residue name and number of ions in pairs, like ['NA',
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
    Run `tleap`, specified by the file and the path.

    Parameters
    ----------
    path : str
        Directory of the `tleap` input file
    file_name : str
        Name of the `tleap` input file

    Returns
    -------
    output : list
        Line-by-line results of running the `tleap` calculation

    """

    utils.check_for_leap_log(path=path)

    p = sp.Popen(['tleap', '-s ', '-f ', file_name], stdout=sp.PIPE, bufsize=1, universal_newlines=True, cwd=path)
    output = []
    # Wait until process terminates...
    while p.poll() is None:
        line = p.communicate()[0]
        output.append(line)
    if p.poll() is None:
        p.kill()
    grep_leap_log(path=path)
    return output

def grep_leap_log(path='./'):
    """
    Check for a few keywords in the `tleap` output.
    """
    with open(path + 'leap.log', r) as file:
        for line in file.readlines():
            if re.search('ERROR|WARNING|Warning|duplicate', line):
                log.warning('It appears there was a problem with solvation: check `leap.log`...')

def basic_tleap(input_file='tleap.in', input_path='./', output_prefix='solvate', output_path=None, pdb_file=None):
    """
    Run `tleap` with a user supplied input file and optionally override the `loadpdb` section. This is usefully for
    quickly solvating a structure without iteratively finding an exact number of water.

    Parameters
    ----------
    input_file : str
        `tleap` input file (this needs to exist already)
    input_path : str
        Directory of the input file
    output_prefix : str
        The name of the re-written input `tleap` input file
    output_path : bool or str
        Directory of the output file (if not specified, will be the same as the input directory)
    pdb_file : bool or str
        Name of `pdb` file to load in `tleap` script

    """

    if output_path is None:
        output_path = input_path
    log.debug('Reading {}/{}, writing {}/{}.in, and executing ...'.format(input_path, input_file, output_path,
                                                                          output_prefix))

    lines = read_tleap_lines(pdb_file=pdb_file, path=input_path, file_name=input_file)
    options = default_tleap_options()
    options['path'] = output_path
    options['output_prefix'] = output_prefix
    write_tleapin(lines, options)
    run_tleap(path=output_path, file_name=output_prefix + '.in')


def count_residues(path='./', file_name='tleap.in'):
    """
    Run and parse `tleap` output and return a dictionary of residues in the structure.
    Parameters
    ----------
    path : str
        Directory of `tleap` input file
    file_name : str
        Name of `tleap` input file

    Returns
    -------
    residues : dict
        Dictionary of added residues and their number
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
    """
    Run and parse `tleap` output and return the volume of the structure.
    Parameters
    ----------
    path : str
        Directory of `tleap` input file
    file_name : str
        Name of `tleap` input file

    Returns
    -------
    volume : float
        The volume of the structure in cubic angstroms

    """
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
    Run and parse `tleap` output and return the a list of water residues.
    Parameters
    ----------
    path : str
        Directory of `tleap` input file
    file_name : str
        Name of `tleap` input file

    Returns
    -------
    water_resiudes : list
        A list of the water residues in the structure
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
    """
    Determine the target number of waters by parsing the `buffer_target` option.

    Parameters
    ----------
    lines : list
        Boilerplate `tleap` input used to get a quick check for a given buffer size in angstroms
    options : dict
        Dictionary containing `tleap` options, set up by `default_tleap_options()` and customized
    buffer_target : str or int
        If `str`, treated as a distance, if `int`, treated as the number of waters

    Returns
    -------
    buffer_target : int
        The desired number of waters for the solvation

    """
    # If buffer_water ends with 'A', meaning it is a buffer distance...
    if isinstance(buffer_target, str) and buffer_target[-1] == 'A':
        # Let's get a rough value of the number of waters if the buffer target is given as a string.
        # This could fail if there is a space in `buffer_target`...
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
    """
    Determine whether additional ions (instead of or in addition to neutralization) are requested...

    Parameters
    ----------
    add_ions : list
        The raw input list passed to the `solvate` function
    options : dict
        Dictionary containing `tleap` options, set up by `default_tleap_options()` or customized
    buffer_target : int
        The desired number of waters (used to calculate molality)
    Returns
    -------
    add_ion_residues : list
        A processed list of ions and their amounts that can be passed to `tleap`

    """
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
            number_to_add = int(np.ceil(number_of_atoms * liters))
            add_ion_residues.append(number_to_add)
        else:
            raise Exception('Unanticipated error calculating how many ions to add.')
    return add_ion_residues


def adjust_buffer_value(number_of_waters, target_number_of_waters, buffer_values, exponent):
    """
    Determine whether to increase or decrease the buffer thickness to match a desired number of waters.

    Parameters
    ----------
    number_of_waters : list
        A list of the number of waters in the solvated structure (see `buffer_values`)
    target_number_of_waters : int
        The desired amount of water
    buffer_values : list
        A list of the number of buffer sizes tried (see `number_of_waters`)
    exponent : float
        A sliding dial that helps narrow dow increases and decreases in buffer thickness

    Returns
    -------
    buffer_value : float
         A new buffer size to try

    """

    # If the last two rounds of solvation have too many waters, make the buffer smaller...
    if number_of_waters[-2] > target_number_of_waters and number_of_waters[-1] > target_number_of_waters:
        # log.debug('Adjustment loop 1')
        return buffer_values[-1] + -1 * (10 ** exponent), exponent

    # If the last two rounds of solvation had too few waters, make the buffer bigger...
    elif number_of_waters[-2] < target_number_of_waters and number_of_waters[-1] < target_number_of_waters:
        # log.debug('Adjustment loop 2')
        return buffer_values[-1] + 1 * (10 ** exponent), exponent

    # If the number of waters was greater than the target and is now less than the target, make the buffer a bit
    # bigger, by an increasingly smaller amount...
    elif number_of_waters[-2] > target_number_of_waters and number_of_waters[-1] < target_number_of_waters:
        # log.debug('Adjustment loop 3')
        exponent -= 1
        return buffer_values[-1] + 5 * (10 ** exponent), exponent

    # If the number of waters was less than the target and is now greater than the target, make the buffer a bit
    # smaller, by an increasingly bigger amount...
    elif number_of_waters[-2] < target_number_of_waters and number_of_waters[-1] > target_number_of_waters:
        # log.debug('Adjustment loop 4')
        exponent -= 1
        return buffer_values[-1] + -5 * (10 ** exponent), exponent
    else:
        raise Exception("The buffer_values search died due to an unanticipated set of variable values")


def remove_waters_manually(lines, number_of_waters, target_number_of_waters, options):
    """
    Remove a few water molecules manually with `tleap` to exactly match a desired number of waters.

    Parameters
    ----------
    lines : str
        Boilerplate `tleap` input read from a template file
    number_of_waters : list
        A list of the number of waters in the solvated structure
    target_number_of_waters : int
        The desired number of waters in the solvated structure
    options : dict
        Dictionary containing `tleap` options, set up by `default_tleap_options()` or customized
    """

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


def solvate(tleap_file, pdb_file=None,
            pbc_type='cubic',
            buffer_target='12.0A',
            water_box='TIP3PBOX',
            neutralize=True,
            counter_cation='Na+',
            counter_anion='Cl-',
            add_ions=None,
            output_prefix='solvate',
            path='./'):
    """
    Solvate a structure with `tleap`, from a template file, to an exact number of waters or buffer size.

    Parameters
    ----------
    tleap_file : str
        Name of template file to read boilerplate `tleap` input (e.g., `source leaprc`... lines)
    pdb_file : str
        Name of `pdb` file to load with `loadpdb`... if not specified in the template file
    pbc_type : str
        Type of solvation (cubic, rectangular, octahedral, or None)
    buffer_target : int or str
        The desired level of water
    water_box : str
        The "type" of water (see AMBER manual for acceptable values)
    neutralize : bool
        Whether to neutralize the system
    counter_cation : str
        If `neutralize`: the positive ion used
    counter_anion : str
        If `neuralize`: the negative ion used
    add_ions : list
        A list of additional ions that can be added to the system, specified by their amount, molarity, or molality
    output_prefix : str
        Name of the written files
    path : str
        Directory of the `tleap_file`

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

        log.debug('Cycle %02d\t %d %d (%d)' % (cycle, options['buffer_value'], waters, target_number_of_waters))

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
            options['buffer_value'], exponent = adjust_buffer_value(number_of_waters, target_number_of_waters,
                                                                    buffer_values,
                                                                    exponent)
            # Now that we're close, let's re-evaluate how many ions to add, in case the volume has changed a lot.
            # (This could be slow and run less frequently...)
            if add_ions and cycle % 10 == 0:
                options['add_ion_residues'] = set_additional_ions(add_ions, options, target_number_of_waters)
            cycle += 1

    if cycle >= max_cycles and waters > target_number_of_waters:
        remove_waters_manually(lines, number_of_waters, target_number_of_waters, options)

    if cycle >= max_cycles and waters < target_number_of_waters:
        raise Exception("Automatic adjustment of the buffer value resulted in fewer waters \
            added than targeted by `buffer_water`. Try increasing the tolerance in the above loop")
    else:
        raise Exception("Automatic adjustment of the buffer value was unable to converge on \
            a solution with sufficient tolerance")
