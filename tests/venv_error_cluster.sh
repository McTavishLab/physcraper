virtualenv -p /usr/bin/python venv-test
source venv-test/bin/activate
git clone https://github.com/McTavishLab/physcraper.git
cd physcraper
pip install --upgrade setupools
python setup.py develop 