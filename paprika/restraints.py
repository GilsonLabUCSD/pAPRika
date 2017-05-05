

# Page 379 of Amber 7 manual has mask specifications
class DAT_restraint(object):
    """
    Distance or angle or torsion restraints.
    """

    def __init__(self):
        self.mask1 = None
        self.mask2 = None
        self.mask3 = None
        self.mask4 = None

        attach =   {'target_intial': None,
                    'target_final':  None,
                    'force_intial':  None,
                    'force_final':   None
                   }
        pull =    {'target_intial':  None,
                   'target_final':   None,
                   'force_intial':   None,
                   'force_final':    None
                  }
        release = {'target_intial':  None,
                   'target_final':   None,
                   'force_intial':   None,
                   'force_final':    None
                  }
        self.phase = {'attach':  attach,
                      'pull':    pull,
                      'release': release}

class POS_restraint(object):
    """
    Positional restraints.
    """

    def __init__(self):
        self.mask1 = None
        self.mask2 = None
        self.mask3 = None
        self.mask4 = None

        attach =   {'target_initial':   None,
                    'target_final':     None,
                    'target_increment': None,
                    'force_initial':    None,
                    'force_final':      None,
                    'force_increment':  None
                   }
        pull =    {'target_initial':   None,
                   'target_final':     None,
                   'target_increment': None,
                   'force_initial':    None,
                   'force_final':      None,
                   'force_increment':  None
                  }
        release = {'target_initial':   None,
                   'target_final':     None,
                   'target_increment': None,
                   'force_initial':    None,
                   'force_final':      None,
                   'force_increment':  None
                  }
        self.phase = {'attach':  attach,
                      'pull':    pull,
                      'release': release}

this = DAT_restraint()
this.mask1 = ':BUT'
this.mask2 = ':CB6'
this.pull = {'target_initial'   : 0,  # angstroms (by definition, I guess)
             'target_final'     : 10, # angstroms
             'target_increment' : 1,  # angstroms 
             'force_intial'     : 5,  # kcal per mol
             'force_final'      : 5,  # kcal per mol
             'force_increment'  : 1   # kcal per mol (but shouldn't matter)
            }
# Generate windows using the target distances
# - Align once, link to each window directory
# - Minimize at each target distance
# - Solvate at each target distance