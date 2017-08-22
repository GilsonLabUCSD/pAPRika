import re as re
import os as os
import numpy as np
import subprocess as sp
import pytraj as pt
import logging as log


# Coordinate Translate
def translate():
    return 0

def gb_minimize(
        tleapfile='tleap.in',
        pdbfile=None,
        restraintfile='restraints.in',
        ntrmask=None,
        saveprefix='vac'):
    """
    Minimize system in implicit solvent, initially with non-bonded off, then gradually non-bondeds turned on.
    """

    # Read tleapfile (mostly duplicated from solvate.py)
    tleaplines = []
    with open(tleapfile, 'r') as f:
        for line in f.readlines():
            if re.search('loadpdb', line):
                words = line.rstrip().replace('=', ' ').split()
                if pdbfile is None:
                    pdbfile = words[2]
                unit = words[0]
                tleaplines.append(
                    "{} = loadpdb {}\n".format(unit, pdbfile))
            if not re.search(r"^\s*addions|^\s*addions2|^\s*addionsrand|^\s*desc|"
                             r"^\s*quit|^\s*solvate|loadpdb|^\s*save", line, re.IGNORECASE):
                tleaplines.append(line)
        tleaplines.append("saveamberparm {0} {1}.prmtop {1}.rst7\n".format(unit,saveprefix))
        tleaplines.append("quit\n")

    # Write tleap file and run
    with open('tleap_gb_minimize.in', 'w') as f:
        for line in tleaplines:
            f.write(line)
    sp.call('tleap -s -f tleap_gb_minimize.in >& log.tleap_gb_minimize', shell=True)

    # Write pmemd input file
    if ntrmask:
        ntrstring="\n  ntr = 1,\n  restraint_wt = 10.0,\n  restraintmask = '{}',".format(ntrmask)
    else:
        ntrstring=""

    with open('gb_minimize.in', 'w') as f:
        f.write("""
*********************************
***** gb_minimize.in ************
*********************************
Minimizing in GB.
 &cntrl
  imin = 1,
  ntx = 1, 
  ntpr = 100,
  maxcyc = 1000,
  ncyc = 250,
  ntxo = 1,
  irest = 0,
  ntf = 1,
  ntc = 1,
  ntb = 0,
  igb = 1,
  cut = 999.0,
  nmropt = 1,
  pencut = -1,{}
 /
 &wt type = 'NB',      istep1=0,    istep2=500,  value1 = 0.000, value2=0.000, IINC=50, /
 &wt type = 'NB',      istep1=500, istep2=750,  value1 = 0.000, value2=1.000, IINC=50, /
 &wt type = 'END', /
DISANG={}
LISTOUT=POUT
        """.format(ntrstring,restraintfile))

    # Run pmemd
    sp.call("pmemd -O -p {0}.prmtop -c {0}.rst7 -ref {0}.rst7 -i gb_minimize.in -o gb_minimize.out -r gb_minimize.rst7 -inf /dev/null".format(saveprefix), shell=True) # Any problems with /dev/null here?

    gbmin = pt.load('gb_minimize.rst7', top=saveprefix+'.prmtop')

    return gbmin
