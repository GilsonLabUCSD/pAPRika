

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

class POS_restraint():
    """
    Positional restraints.
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