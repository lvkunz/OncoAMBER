python setup.py sdist
rm -rf OncoAMBER.egg-info
mv dist/* .
rm -rf dist/
conda install -c conda-forge OncoAMBER-1.1.3.tar.gz
