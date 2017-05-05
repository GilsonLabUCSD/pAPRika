import os
import sys
import versioneer
from distutils.core import setup

if sys.version_info < (2, 7):
    sys.stderr.write('You must have at least Python 2.7 for pAPRika to work '
                     'correctly, although we have only tested with Python 3.0.\n')
    sys.exit(1)

from distutils.command.clean import clean as Clean

class CleanCommand(Clean):
    """python setup.py clean
    """
    # lightly adapted from scikit-learn package
    description = "Remove build artifacts from the source tree"

    def _clean(self, folder):
        for dirpath, dirnames, filenames in os.walk(folder):
            for filename in filenames:
                if (filename.endswith('.so') or filename.endswith('.pyd')
                        or filename.endswith('.dll')
                        or filename.endswith('.pyc')):
                    os.unlink(os.path.join(dirpath, filename))
            for dirname in dirnames:
                if dirname == '__pycache__':
                    shutil.rmtree(os.path.join(dirpath, dirname))

    def run(self):
        Clean.run(self)
        if os.path.exists('build'):
            shutil.rmtree('build')
        self._clean('./')

# paprika package and all its subpackages
packages = ['paprika']


if __name__ == '__main__':

    cmdclass = dict(clean=CleanCommand)
    cmdclass.update(versioneer.get_cmdclass())
    setup(name='pAPRika',
          version=versioneer.get_version(),
          description='Attach-pull-release free energy calculations',
          author='',
          author_email='',
          url='',
          license='',
          packages=packages,
          cmdclass=cmdclass)
