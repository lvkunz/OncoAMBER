
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import io
import os
import sys
import re
from shutil import rmtree
from distutils.core import setup
from setuptools import find_packages, Command

# Package meta-data
NAME = 'OncoAMBER'
DESCRIPTION = 'Agent-based model of tumor growth and response to radiation therapy'
URL = 'https://github.com/lvkunz/OncoAMBER'
EMAIL = 'lvkunz@mgh.harvard.edu'
AUTHOR = 'Louis Kunz'
REQUIRES_PYTHON = '>=3.8.2'

# Required packages for this module to be executed
REQUIRED = ['numpy', 'pandas', 'scipy']

# Optional packages
EXTRAS = { 'plots' :['matplotlib', 'seaborn'] }

here = os.path.abspath(os.path.dirname(__file__))

# Import the README and use it as the long-description
try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except:
    long_description = DESCRIPTION

# Read the version from __init__.py
version = ''
init_file = os.path.join(os.path.dirname(__file__), 'amber', '__init__.py')
with open(init_file, 'r') as f:
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", f.read(), re.M)
    if version_match:
        version = version_match.group(1)

# Assign the version to the VERSION variable
VERSION = version

# Support setup.py upload
class UploadCommand(Command):
    description = 'Build and publish the package'
    user_options = []

    @staticmethod
    def status(s):
        """Prints things in bold."""
        print('\033[1m{0}\033[0m'.format(s))

    def initialize_options(self):
        pass

    def finalize_options(self):
        pass

    def run(self):
        try:
            self.status('Removing previous builds...')
            rmtree(os.path.join(here, 'dist'))
        except OSError:
            pass

        self.status('Building source and wheel distribution...')
        os.system('{0} setup.py sdist bdist_wheel --universal'.format(sys.executable))

        self.status('Uploading the package to PyPI via Twine...')
        os.system('twine upload dist/*')

        self.status('Pushing git tags...')
        os.system('git tag v{0}'.format(VERSION))
        os.system('git push --tags')

# Executing setup
setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(),
    package_data={'': ['*.txt', '*.csv'],},
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    # include_package_data=True,
    zip_safe=False,
    license='MIT',
    classifiers=[   'License :: OSI Approved :: MIT License',    'Operating System :: MacOS',    'Operating System :: Microsoft :: Windows',    'Programming Language :: Python',    'Programming Language :: Python :: 3',    'Programming Language :: Python :: 3.6',    'Programming Language :: Python :: Implementation :: CPython',    'Programming Language :: Python :: Implementation :: PyPy'],
    cmdclass={
        'upload' : UploadCommand,
    }
)