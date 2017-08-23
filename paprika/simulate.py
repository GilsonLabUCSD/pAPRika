import numpy as np
import subprocess as sp
import logging as log
import os as os
import parmed as pmd

class simulation()
    """
    All simulation methods.  Minimization, equilibration, production.
    """

    def __init__(self):

        self.minimization_settings = {
            'imin':         1,

        } 
