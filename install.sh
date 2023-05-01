python setup.py sdist
rm -rf OncoAMBER.egg-info
mv dist/* .
rm -rf dist/
pip install OncoAMBER-1.1.3.tar.gz
