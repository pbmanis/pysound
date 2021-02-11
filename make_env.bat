Rem @echo off
set ENVNAME=ps_venv
set PSPATH=%cd%
py -3.8 -m venv %ENVNAME%
set ACT=C:\%PSPATH%\%ENVNAME%\Scripts\activate.bat
set DEACT=C:\%PSPATH%\%ENVNAME%\Scripts\deactivate.bat
%ACT%
python --version

python -m pip install --upgrade pip 
Rem  be sure pip is up to date in the new env.
pip3 install wheel  
Rem seems to be missing (note singular)
pip3 install cython
Rem # if requirements.txt is not present, create:
Rem # pip install pipreqs
Rem # pipreqs

Rem See noet about how to get pyaudio in the requirements.txt file
pip3 install -r requirements.txt
pip3 install downloads\\PyAudio-0.2.11-cp38-cp38-win_amd64.whl
CALL %ENVNAME%\Scripts\activate.bat
python setup.py develop

Rem Should always run test afterwards.
python tests/play_test_sounds.py

pause