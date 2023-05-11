#!/bin/bash

python setup.py sdist  # Create the source distribution package
rm -rf OncoAMBER.egg-info  # Remove the existing egg-info directory
mv dist/* .  # Move the distribution package to the current directory
rm -rf dist/  # Remove the dist directory

# Install the package
pip install $(find . -name "OncoAMBER-*.tar.gz" | sort -V | tail -n 1)

# Clean up the installation files
rm -f OncoAMBER-*.tar.gz
