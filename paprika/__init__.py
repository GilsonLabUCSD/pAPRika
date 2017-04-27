"""
The pAPRika package sets up and performs attach-pull-release calculations.
"""

from ._version import get_versions
import argparse as argparse
import logging as logging
import parmed as pmd

__version__ = get_versions()['version']
del get_versions


# This is nice, but we should move away from argparse if we're going
# to call this as a module.
parser = argparse.ArgumentParser()
parser.add_argument('-d',
                    '--debug',
                    help='Really increase output verbosity',
                    action='store_const',
                    dest='loglevel',
                    const=logging.DEBUG,
                    default=logging.WARNING)

parser.add_argument('-v', 
                    '--verbose',
                    help='Increase output verbosity',
                    action='store_const',
                    dest='loglevel',
                    const=logging.INFO)
args = parser.parse_args()
# logging.basicConfig(level=args.loglevel)
# logging.basicConfig(format='%(asctime)s %(message)s', datefmt='%Y-%m-%d %I:%M:%S %p',
#                    level=args.loglevel)
