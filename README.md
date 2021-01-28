pysounds
========

A set of modules that provide auditory stimulus generation, including hardware control.

It includes:

sound.py
--------
A collection of routines that generate sound waveforms. Uses DMR.py for dynamic ripple, 

pystim.py
---------
A module that interacts with some hardware (Rp2.1, PA5, RZ5D, nidaq, and soundcards) to play
generated sounds.

pysounds.py
-----------
An early verion module that only knows about the RP2.1, PA5, and nidaq. 

calreader.py
------------
Read and display calibration files generated in matlab ABR program.

stimController
--------------

This module suite provides the following in conjunction with TDT hardware:

1. Generation of waveforms for acoustic stimulation, including tones, clicks, noises, bandpass and notched noise, modulated tones and noises, RSS, DMR, and comodulated masking stimuli (not all are currently implemented).

2. The suite uses a graphical window to show the waveforms and long-term spectra, as well as an option to display the frequency-time spectra for the stimuli. 

3. The suite provides for the presentation of stimuli, including control of their timing, repetitions, etc. 

4. The suite interfaces with system soundcards (for testing), National Instruments DAC card(s)*, and Tucker-Davis Technologies RP2.1, RZ5D DSP processors and PA5 attenuators.

5. For recording responses to acoustic stimulation, the suite expects to connect to an RZ5D with a specific circuit that is loaded. In this setup, the RZ5D (using appropriately configured OpenWorkbench) generates trigger pulses for the stimuli, and records multiple channels from the preamplifier. The current suite in turn controls the RZ5D, starting the acquisition, and following the storage of data in the data tank. The goal is that most of the experimental interaction with the user takes place with the current suite, rather than through the RZ5D/OpenWorkbench/OpenScope windows, although these need to be monitored. Communication with the RZ5D is limited to setting up the standard CoreSweepControl timing parameters, and initiating mode changes (Record, Preview, Standby, Idle). 


This sofware does not provide:

1. Data analysis.

2. Complex stimulus control.



* Tested with NI6371 only.