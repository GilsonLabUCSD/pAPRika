import logging as log
import subprocess as sp
from collections import OrderedDict
import os
import time

class Simulation(object):
    """
    AMBER simulation class.
    """

    def __init__(self):

        ### Setup simulation directory and files
        self.path = '.' # Assume everything will be created/executed in this path
        self.executable = 'sander'
        self.CUDA_VISIBLE_DEVICES = None
        self.phase = None   
        self.window = None
        self.topology = 'prmtop'
        self.restraint_file = 'restraints.in'
        self.title = 'PBC MD Simulation'
        self.converged = False

        ### File names
        self._prefix = 'md'
        self.input = self._prefix + '.in'
        self.inpcrd = self._prefix + '.inpcrd'
        self.ref = self._prefix + '.inpcrd'
        self.output = self._prefix + '.out'
        self.restart = self._prefix + '.rst7'
        self.mdinfo = self._prefix + '.mdinfo'
        self.mdcrd = self._prefix + '.nc'
        self.mden = self._prefix + '.mden'

        ### Input file cntrl settings (Default = NTP)
        self.cntrl = OrderedDict()
        self.cntrl['imin'] = 0
        self.cntrl['ntx'] = 1
        self.cntrl['irest'] = 0
        self.cntrl['maxcyc'] = 0
        self.cntrl['ncyc'] = 0
        self.cntrl['dt'] = 0.002
        self.cntrl['nstlim'] = 5000
        self.cntrl['ntpr'] = 500
        self.cntrl['ntwe'] = 500
        self.cntrl['ntwr'] = 5000
        self.cntrl['ntwx'] = 500
        self.cntrl['ntxo'] = 1
        self.cntrl['ioutfm'] = 1
        self.cntrl['ntf'] = 2
        self.cntrl['ntc'] = 2
        self.cntrl['cut'] = 8.0
        self.cntrl['igb'] = 0
        self.cntrl['tempi'] = 298.15
        self.cntrl['temp0'] = 298.15
        self.cntrl['ntt'] = 3
        self.cntrl['gamma_ln'] = 1.0
        self.cntrl['ig'] = -1
        self.cntrl['ntp'] = 1
        self.cntrl['barostat'] = 2
        self.cntrl['ntr'] = None
        self.cntrl['restraint_wt'] = None
        self.cntrl['restraintmask'] = None
        self.cntrl['nmropt'] = 1
        self.cntrl['pencut'] = -1

        # Other input file sections
        self.ewald = None
        self.other_namelist = None # Could add other namelists as dicts
        self.wt = None    # or []
        self.group = None    # or []


    ### Refresh file names if prefix changes
    @property
    def prefix(self):
        return self._prefix
    @prefix.setter
    def prefix(self, new_prefix):
        self._prefix = new_prefix
        self.input = new_prefix + '.in'
        self.inpcrd = new_prefix + '.inpcrd'
        self.ref = new_prefix + '.inpcrd'
        self.output = new_prefix + '.out'
        self.restart = new_prefix + '.rst7'
        self.mdinfo = new_prefix + '.mdinfo'
        self.mdcrd = new_prefix + '.nc'
        self.mden = new_prefix + '.mden'

    def _config_min(self):
        """
        Configure input settings for minimization.
        """
        self.cntrl['imin'] = 1
        self.cntrl['ntx'] = 1
        self.cntrl['irest'] = 0
        self.cntrl['maxcyc'] = 5000
        self.cntrl['ncyc'] = 1000
        self.cntrl['dt'] = 0.0
        self.cntrl['nstlim'] = 0
        self.cntrl['ntpr'] = 100
        self.cntrl['ntwr'] = 5000
        self.cntrl['ntwx'] = 0
        self.cntrl['ntwe'] = 0
        self.cntrl['ntxo'] = 1
        self.cntrl['ntf'] = 1
        self.cntrl['ntc'] = 1
        self.cntrl['ntt'] = 0
        self.cntrl['gamma_ln'] = 0.0
        self.cntrl['ig'] = 0
        self.cntrl['ntp'] = 0
        self.cntrl['barostat'] = 0
        self.mdcrd = None
        self.mden = None


    def config_pbc_min(self):
        """
        Configure input settings to minimization in periodic boundary conditions.
        """
        self._config_min()
        self.title = 'PBC Minimization'
        self.cntrl['cut'] = 8.0
        self.cntrl['igb'] = 0


    def config_gb_min(self):
        """
        Configure input settings to minimization in continuum solvent.
        """

        self._config_min()
        self.title = 'GB Minimization'
        self.cntrl['cut'] = 999.0
        self.cntrl['igb'] = 1


    def _config_md(self):
        """
        Configure input settings for MD.
        """
        self.cntrl['imin'] = 0
        self.cntrl['ntx'] = 1
        self.cntrl['irest'] = 0
        self.cntrl['maxcyc'] = 0
        self.cntrl['ncyc'] = 0
        self.cntrl['dt'] = 0.002
        self.cntrl['nstlim'] = 5000
        self.cntrl['ntpr'] = 500
        self.cntrl['ntwe'] = 500
        self.cntrl['ntwr'] = 5000
        self.cntrl['ntwx'] = 500
        self.cntrl['ntxo'] = 1
        self.cntrl['ioutfm'] = 1
        self.cntrl['ntf'] = 2
        self.cntrl['ntc'] = 2
        self.cntrl['ntt'] = 3
        self.cntrl['gamma_ln'] = 1.0
        self.cntrl['ig'] = -1

    def config_gb_md(self):
        """
        Configure input settings for MD in default GB.
        """

        self._config_md()
        self.title = 'GB MD Simulation'
        self.cntrl['cut'] = 999.0
        self.cntrl['igb'] = 1
        self.cntrl['ntp'] = 0
        self.cntrl['barostat'] = 0


    def config_pbc_md(self):
        """
        Configure input settings to default NTP.
        """

        self._config_md()
        self.title = 'PBC MD Simulation'
        self.cntrl['cut'] = 8.0
        self.cntrl['igb'] = 0
        self.cntrl['ntp'] = 1
        self.cntrl['barostat'] = 2

    def _write_dict_to_mdin(self, f, dictionary):
        for key, val in dictionary.items():
            if val is not None:
                f.write("  {:15s} {:s},\n".format(key + ' =', str(val)))
        f.write(" /\n")


    def _amber_write_input_file(self):
        log.debug('Writing {}'.format(self.input))
        with open(self.path+'/'+self.input, 'w') as f:
            f.write("{}\n".format(self.title))

            f.write(" &cntrl\n")
            self._write_dict_to_mdin(f, self.cntrl)

            if self.ewald is not None:
                f.write(" &ewald\n")
                self._write_dict_to_mdin(f, self.ewald)

            if self.cntrl['nmropt'] == 1:
                if self.wt is not None:
                    for line in self.wt:
                        f.write(" " + line + "\n")
                f.write(" &wt type = 'END', /\n")
                if self.restraint_file is not None:
                    f.write("DISANG = {}\n".format(self.restraint_file))
                    f.write("LISTOUT = POUT\n\n")
            if self.group is not None:
                f.write("{:s}".format(self.group))

    
    def run(self, soft_minimize=False, overwrite=False):
        """
        Minimize the system.

        If soft=True, slowly turn on non-bonded interactions during minimization
        so that the restraints get enforced first.

        """

        if overwrite or not self.has_timings():

            # These settings hardcoded at the moment ... possibly expose for editing in the future
            if soft_minimize:
                # Set a burn in value that is 25% of the way between ncyc and maxcyc
                ncyc = self.cntrl['ncyc']
                maxcyc = self.cntrl['maxcyc']
                burn_in = int(float(ncyc) + 0.20*(float(maxcyc) - float(ncyc)))
                # If the burn_in value is nuts, then just set it to zero
                if burn_in < 0 or burn_in >= maxcyc:
                    burn_in = 0
                # Set an end_soft value that is 75% of way between ncyc and maxcyc
                end_soft = int(float(ncyc) + 0.60*(float(maxcyc) - float(ncyc)))
                self.wt = [
                    "&wt type = 'NB', istep1=0, istep2={:.0f}, value1 = 0.0, value2=0.0, IINC=50, /".format(burn_in),
                    "&wt type = 'NB', istep1={:.0f}, istep2={:.0f}, value1 = 0.0, value2=1.0, IINC=50, /".format(burn_in,end_soft)]
    
            #_amber_write_input_file(self.path+'/'+self.input, self.min, title='GB Minimization.')
            self._amber_write_input_file()
    
            if self.cntrl['imin'] == 1:
                log.info('Running Minimization at {}'.format(self.path))
            else:
                log.info('Running MD at {}'.format(self.path))
    
            # Create executable list for subprocess
            exec_list = self.executable.split() + ['-O', '-p', self.topology]
            if self.ref is not None:
                exec_list += ['-ref', self.ref]
            exec_list += ['-c', self.inpcrd, '-i', self.input, '-o', self.output, '-r', self.restart]
            if self.mdcrd is not None:
                exec_list += ['-x', self.mdcrd]
            if self.mdinfo is not None:
                exec_list += ['-inf', self.mdinfo]
            if self.mden is not None:
                exec_list += ['-e', self.mden]

            log.debug('Exec line: '+' '.join(exec_list))

            # Execute
            if self.CUDA_VISIBLE_DEVICES:
                amber_output = sp.Popen(exec_list, cwd=self.path, stdout=sp.PIPE, stderr=sp.PIPE,
                                        env=dict(os.environ, CUDA_VISIBLE_DEVICES=str(self.CUDA_VISIBLE_DEVICES)))
            else:
                amber_output = sp.Popen(exec_list, cwd=self.path, stdout=sp.PIPE, stderr=sp.PIPE)

            amber_output = amber_output.stdout.read().splitlines()

            # Report any stdout/stderr which are output from execution
            if amber_output:
                log.info('STDOUT/STDERR received from AMBER execution')
                for line in amber_output:
                    log.info(line)

            # Check completion status
            if self.cntrl['imin'] == 1 and self.has_timings():
                log.info('Minimization completed...')
            elif self.has_timings():
                log.info('MD completed ...')
            else:
                log.info('Simulation execution does not appear to have completed')

        else:
            log.info("Completed output detected ... Skipping. Use: run(overwrite=True) to overwrite")

    def has_timings(self, alternate_file=None):
        """
        Check for the string TIMINGS in self.ouput file.

        Parameters
        ----------
        alternate_file : str
            If present, check for TIMINGS in this file rather than self.output. Default: None

        Returns
        -------
        timings : bool
            True if 'TIMINGS' is found in file. False, otherwise. 

        """

        # Assume not completed
        timings = False

        if alternate_file:
            output_file = alternate_file
        else:
            output_file = os.path.join(self.path, self.output)

        if os.path.isfile(output_file):
            with open(output_file, 'r') as f:
                strings = f.read()
                if (' TIMINGS' in strings):
                    timings = True

        return timings


