python setup.py sdist
rm -r OncoAMBER.egg-info
mv dist/* .
rm -r dist/
pip install OncoAMBER-1.1.0.tar.g