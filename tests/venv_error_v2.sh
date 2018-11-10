rm -r docs/source/
rm docs/Makefile
rm make.bat

virtualenv -p /usr/bin/python venv-test2
source venv-test2/bin/activate
python setup.py develop 