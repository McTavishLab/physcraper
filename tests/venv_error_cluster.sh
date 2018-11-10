virtualenv -p /usr/bin/python venv-test
source venv-test/bin/activate
git clone https://github.com/McTavishLab/physcraper.git
cd physcraper
pip install --upgrade setuptools
python setup.py develop
pip install -r requirements.txt
python tests/testfilesetup.py
sh tests/run_tests.sh