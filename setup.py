import os
import sys
import re
import codecs

from setuptools import setup, find_packages

if sys.version_info < (2, 7):
    sys.stderr.write(
        'You must have at least Python 2.7 for pAPRika to work '
        'correctly, although we have only tested with Python 3.0.\n')
    sys.exit(1)

# https://packaging.python.org/guides/single-sourcing-package-version/#single-sourcing-the-version
here = os.path.abspath(os.path.dirname(__file__))


def read(*parts):
    with codecs.open(os.path.join(here, *parts), 'r') as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]",
                              version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


if __name__ == '__main__':

    setup(
        name='pAPRika',
        version=find_version("paprika", "__init__.py"),
        description='Attach-pull-release free energy calculations',
        author='',
        author_email='',
        url='',
        license='',
        packages=['paprika'])
