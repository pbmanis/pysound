ENVNAME="ps_venv"

rm -rf $ENVNAME
py -3.10 -m venv $ENVNAME
source $ENVNAME/Scripts/activate

python -m pip install --upgrade pip  # be sure pip is up to date in the new env.
pip3 install wheel  # seems to be missing (note singular)
pip3 install cython>=3.0

# # if requirements.txt is not present, create:
# # pip install pipreqs
# # pipreqs
#
# #Then:
#
pip3 install -r requirements.txt
source $ENVNAME/Scripts/activate
python Scripts/pywin32_postinstall.py -install

python --version
python setup.py develop

