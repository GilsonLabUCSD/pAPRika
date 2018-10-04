import os
import sys

from setuptools import setup
from paprika.version import find_version

if sys.version_info < (2, 7):
    sys.stderr.write('You must have at least Python 2.7 for pAPRika to work '
                     'correctly, although we have only tested with Python 3.0.\n')
    sys.exit(1)

# https://packaging.python.org/guides/single-sourcing-package-version/#single-sourcing-the-version
here = os.path.abspath(os.path.dirname(__file__))

if __name__ == '__main__':

    setup(
        name='pAPRika',
        version=find_version(),
        description='Attach-pull-release free energy calculations',
        author='',
        author_email='',
        url='',
        license='',
        packages=['paprika'])
