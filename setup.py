import os
import sys
import re
from distutils.core import setup
import subprocess as sp

if sys.version_info < (2, 7):
    sys.stderr.write(
        'You must have at least Python 2.7 for pAPRika to work '
        'correctly, although we have only tested with Python 3.0.\n')
    sys.exit(1)

version = '0.0.2'
try:
    # use git to find current version
    git_describe = sp.check_output(["git", "describe", "--always"]).strip()
    version = re.sub('-g[0-9a-f]*$', '', git_describe)
except:
    pass

if __name__ == '__main__':

    setup(
        name='pAPRika',
        version=version,
        description='Attach-pull-release free energy calculations',
        author='',
        author_email='',
        url='',
        license='',
        packages='paprika')
