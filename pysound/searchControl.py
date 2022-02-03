"""
Stimulus controller
Generates waveforms, controls stimulus presentation
Reiles on sound.py for stimulus waveform generation
Relies on pystim.py for hardware interactions.

Operates in two modes for output intensity:

Atten mode: the levels refer to attenuation in Db, no correction for stimlus
SPL mode : the levels are corrected by the system calibration on a per frequency basis
Currently, only the atten mode is supported. 

"""

######################- Stimulus repetition code- need to decide where this goes...  Big Idea- load the DAQ once and retrigger
# stimulus_duration = isi * reps  # len(wavel)*samplefreq + postduration
# pts_per_rep = int(float(isi) * samplefreq)
# if wavel.shape[0] < pts_per_rep:
#     wavel = np.concatenate(
#         (wavel, np.zeros(pts_per_rep - wavel.shape[0])), axis=0
#     )
# wavel = np.tile(wavel, reps)
# if waver is not None:
#     if waver.shape[0] < pts_per_rep:
#         waver = np.concatenate(
#             (waver, np.zeros(pts_per_rep - waver.shape[0])), axis=0
#         )
#     waver = np.tile(waver, reps)

###########################---End Stimulus repetition code

import sys
import os
from dataclasses import dataclass, field
import datetime
import importlib
import numpy as np
import scipy.signal
import time
import pickle
import scipy.io.wavfile as wav
from collections import OrderedDict

# from backports import configparser
import pyqtgraph as pg
from PyQt5 import QtGui, QtCore
from PyQt5.QtCore import pyqtSlot, pyqtSignal, QThreadPool, QTimer, QRunnable
from pyqtgraph.parametertree import Parameter, ParameterTree
from pysound import pystim
from pysound import Utility
from pysound import sound
import pprint

# import TDTTankInterface as TDT
import tdt  # necessary for interacting with the tanks in Synapse

pp = pprint.PrettyPrinter(indent=4)

@dataclass
class ControllerState:
    """
    All flags related to the controller state belong here
    Includes:
        Current status of RZ5D (Idle, Preview, Running)
        Status of NI DAC output (running, not running)
        # of stimuli since last NI stimulus load ("repetitions")
        Runtime since last NI load
        Runtime since last Record start
    """
    RZ5DStatus: str = "Idle"
    NIDAQStatus: str = "stopped"
    N_stimuli: int = 0
    Stimulus_runtime: float = 0.
    Recording_runtime: float = 0.

class Controller(object):
    """
    Controller provides control over the stimlus hardware and data acquisition
    """
    def __init__(self, ptreedata, plots, img, maingui):
        self.ptreedata = ptreedata
        self.plots = plots  # access to plotting area
        self.img = img
        self.status = ControllerState()

        self.maingui = maingui
        # self.ProtocolNumber = 0
        self.TDTinfo = tdt.SynapseAPI()
        self.getAllParameters(ptreedata)  # read the current stimlus parameters

        self.PS = pystim.PyStim(
            required_hardware=["PA5", "NIDAQ", "RZ5D"], controller=self,
        )  # , 'RP21'])

        # Special for clickable map
        # we don't save the data so we don't directly program these - they change with the point clicked in the map

    # --------TFR---- we don't need these
    # self.attn = 35
    # self.tone_frequency = 4.0 # khz

    def getAllParameters(self, params):
        """
        Get all of the local parameters from the parameter tree
        and set into "Cpars"

        Parameters
        ----------
        ptree : ParameterTree object

        Returns
        -------
        Nothing
        """
        # fill the Parameter dictionary from the parametertree
        self.CPars = OrderedDict()
        self.CPars["IntensityMode"] = "attenuation"
        self.CPars["Voltage_Scales"] = {
            "Tone_V": 10.0,
            "maxTone_dB": {"MF1": 110, "EC1": 83.9},
            "Click_V": 5.0,
            "maxClick_dB": {"MF1": 108.5, "EC1": 79.5},
            "Noise_V": 2.5,
            "maxNoise_dB": {"MF1": 0, "EC1": 0},
        }  # we don't actualy know... but also need to clip

        for ch in self.ptreedata.childs:
            self.CPars[ch.name()] = {}
            for par in ch.childs:
                # print(' name: %s ' % par.name()),
                if par.type() == "int":
                    self.CPars[ch.name()][par.name()] = int(par.value())
                elif par.type() == "float":
                    self.CPars[ch.name()][par.name()] = float(par.value())
                elif par.type() == "list":
                    self.CPars[ch.name()][par.name()] = str(par.value())
                elif par.type() == "str":
                    self.CPars[ch.name()][par.name()] = str(par.value())

    def change(self, param, changes):
        """
        Respond to changes in the parametertree and update class variables
        
        Parameters
        ----------
        param : parameter list
        
        changes : changes as returned from the parameter tree object
        
        Returns
        -------
        Nothing
        
        """
        for param, change, data in changes:
            path = self.ptreedata.childPath(param)
            self.CPars[path[0]][path[1]] = data  # 2 levels only...
        # When a change occurs, we should also update the stimulus if it is running.
        self.stop_stim()
        self.prepare_stimulus()
        print(self.PS.get_RZ5D_Mode())
        if self.PS.get_RZ5D_Mode() in ['Preview']:
            self.start_run_preview()
        elif self.PS.get_RZ5D_Mode() in ['Record']:
            self.start_run_record()

        # self.showParameters()

    def showParameters(self):
        for k in list(self.CPars.keys()):
            print(("Group: %s" % k))
            if isinstance(self.CPars[k], dict):
                for d in list(self.CPars[k].keys()):
                    print(("   %s = %s" % (d, str(self.CPars[k][d]))))

    def show_stimulus_count(self):
        self.status.N_stimuli = self.PS.State.stimulus_count
        self.maingui.stimulus_counter.setText(f"# Stim: {self.status.N_stimuli:4d}")

    def start_run_preview(self):
        """
        Initialize variables for the start of a run
        then use next_stimulus to compute waveform and start timers
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Nothing
        """
        self.clearErrMsg()

        self.runtime = 0
        self.status.N_stimuli = 0
        self.startTime = datetime.datetime.now()
        self.maingui.label_status.setText("Preview")
        self.maingui.label_status.setStyleSheet("font:bold 14px; color: green")

        self.prepare_stimulus()  # reset the data arrays and calculate the next stimulus
        self.PS.play_sound(
                wavel=self.wave,
                waver=self.wave,
                samplefreq=self.PS.Stimulus.out_sampleFreq,
                isi=self.CPars["Stimulus"]["Interstimulus Interval"],
                reps=self.CPars["Stimulus"]["Max Repetitions"],
                attns=self.convert_spl_attn(self.CPars["Stimulus"]["Attenuator"]),
                storedata=False,
            )

    def start_run_record(self):
        """
        Initialize variables for the start of a run
        then use next_stimulus to compute waveform and start timers
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Nothing
        """
        self.clearErrMsg()

        self.runtime = 0
        self.status.N_stimuli = 0
        self.startTime = datetime.datetime.now()
        self.maingui.label_status.setText("Recording")
        self.maingui.label_status.setStyleSheet("font:bold 14px; color: red")

        self.prepare_stimulus()  # reset the data arrays and calculate the next stimulus
        self.PS.play_sound(
                wavel=self.wave,
                waver=self.wave,
                samplefreq=self.PS.Stimulus.out_sampleFreq,
                isi=self.CPars["Stimulus"]["Interstimulus Interval"],
                reps=self.CPars["Stimulus"]["Max Repetitions"],
                attns=self.convert_spl_attn(self.CPars["Stimulus"]["Attenuator"]),
                storedata=True,
            )


    def stop_stim(self):
        """
        Turn off the stimulus; but leave the acquisition running
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Nothing
        """
        self.PS.stop_nidaq()
        self.maingui.label_status.setText("NI Paused")
        self.maingui.label_status.setStyleSheet("font:bold 14px; color: cyan")


    def stop_run(self):
        """
        Stop all stimuli and acquisition.
        Write data
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Nothing
        """
        self.PS.stop_recording()
        self.maingui.label_status.setText("Stopped")
        self.maingui.label_status.setStyleSheet("font:bold 14px; color: black")

        return

    def quit(self):
        """
        Provide a clean exit, by stopping all stimuli and disengaging hardware
        """
        self.PS.stop_recording()
        self.PS.HwOff()
        exit(0)

    def convert_spl_attn(self, spl):
        if self.CPars["IntensityMode"] == "attenuation":
            return spl  # use as attenuation directly
        elif self.CPars["IntensityMode"] == "spl":
            return [
                100.0 - spl,
                100.0 - spl,
            ]  # just rough, within 15dB, need to clean up
        else:
            raise ValueError("Intensity mode must be attenuation or spl")

    def map_voltage(self, protocol, wave, clip=True):
        """
        Provide scaling of the stimulus voltage based on the stimuli. 
        General rule: use the maximal voltage as possible while avoiding
        the possiblity of clipping. Let the attenuators do the rest of the work.
        
        Parameters
        ----------
        protocol : str (no default)
            The name of the protocol
        
        wave : numpy array (no default)
            The 1D data array that will be scaled
        
        clip : boolean (default: True)
            Set true to clip the waveform at +/- 10 V
        
        """
        knownprotocols = ["Noise Search", "Tone Search", "Click Search"]
        if protocol not in knownprotocols:
            raise ValueError(
                "Protocol not in list we can map for scaling the voltage in map_voltage"
            )
        if protocol in ["Click Search"]:
            A = self.CPars["Voltage_Scales"]["Click_V"]
        else:
            A = 10.0

        waves = wave * A
        if clip:
            waves[waves > 10.0] = 10.0
            waves[waves < -10.0] = -10.0
        return waves

    def prepare_stimulus(self, freq=None, level=None):
        """
        Compute the new stimulus
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Nothing
        """
        stim = self.CPars["Stimulus"]["Protocol"]
        #        level = None  # level is dbspl normally for models, but set to None for TDT (5V tone reference)

        wave = None
        # self.stim_vary = None
        # self.total_trials = 1000
        if stim in ["Click Search"]:
            wave = sound.ClickTrain(
                rate=self.PS.Stimulus.out_sampleFreq,  # sample frequency,
                duration=self.CPars["Stimulus"]["Duration"] + self.CPars["Stimulus"]["Delay"],
                dbspl=level,
                click_duration=self.CPars["Stimulus"]["Click Duration"]*1e-6,
                click_starts=np.arange(0.1, 1),
            )

        elif stim in ["Tone Search"]:
            if freq is None:
                freq = self.CPars["Stimulus"]["Tone Frequency"] * 1000.0
            wave = sound.TonePip(
                rate=self.PS.Stimulus.out_sampleFreq,  # sample frequency,
                duration = self.CPars["Stimulus"]["Duration"] + self.CPars["Stimulus"]["Delay"],
                f0=freq,
                dbspl=level,
                pip_duration=self.CPars["Stimulus"]["Duration"],
                # - self.CPars["Stimulus"]["Interstimulus Interval"],
                pip_starts=[self.CPars["Stimulus"]["Delay"]],
                ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000,
            )

        elif stim in ["Noise Search"]:
            seed = 32767
            wave = sound.NoisePip(
                rate=self.PS.Stimulus.out_sampleFreq,  # sample frequency,
                duration=self.CPars["Stimulus"]["Duration"] + self.CPars["Stimulus"]["Delay"],
                dbspl=level,
                pip_duration=self.CPars["Stimulus"]["Duration"],
                pip_starts=[self.CPars["Stimulus"]["Delay"]],
                ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000,
                seed=seed,
            )
            # wave = sound.NoisePip(rate=Fs, duration=self.CPars['Stimulus']['Duration'],
            #             f0=self.CPars['Stimulus']['Tone Frequency']*1000., dbspl=level,
            #             pip_duration=self.CPars['Stimulus']['Duration']-self.CPars['Stimulus']['Interstimulus Interval'],
            #             pip_start=np.arange(0.1,1), ramp_duration=self.CPars['Stimulus']['Rise-Fall']/1000,seed=seed)

        if wave is not None:
            self.wavesound = wave
            self.wave = self.map_voltage(
                stim, self.wavesound.sound, clip=True
            )  # force computation and rescale and clip the waveform

    def show_wave(self):
        """
        Plot the waveform in the top graph
        """
        self.clearErrMsg()
        self.prepare_stimulus()  # force computation/setup of stimulus
        self.plots["Wave"].clear()
        # self.plots['Wave'].plot(self.wavesound.time, self.wave)
        self.plots["Wave"].plot(self.wavesound.time, self.wave)

    def show_spectrogram(self):
        """
        Plot the spectrum in the middle graph
        
        If spectimage is checked in the main gui, also plot the spectrogram
        in a matplotlib window (must be closed to continue; is blocking)
        Note: we have commented the matplotlib plot out, just use pyqtgraph
        """
        self.clearErrMsg()
        self.prepare_stimulus()
        # show the long term spectrum.

        f, Pxx_spec = scipy.signal.periodogram(
            self.wave, self.PS.Stimulus.out_sampleFreq
        )  # , window='flattop', nperseg=8192,
        # noverlap=512, scaling='spectrum')
        self.plots["LongTermSpec"].clear()
        self.plots["LongTermSpec"].plot(f[1:], np.sqrt(Pxx_spec)[1:], pen=pg.mkPen("y"))
        # self.plots['LongTermSpec'].setLogMode(x=True, y=False)

        # f, Pxx_spec = scipy.signal.periodogram(self.wave, Fs) #, window='flattop', nperseg=8192,
        #                # noverlap=512, scaling='spectrum')
        # self.plots['LongTermSpec'].clear()
        # self.plots['LongTermSpec'].plot(f[1:], np.sqrt(Pxx_spec)[1:], pen=pg.mkPen('y'))

    def clearErrMsg(self):
        """
        Reset the error notificatoin to the standard ready indicator
        """
        self.maingui.permStatusMessage.setText('<b><font color="#00FF00">Ready</b>')


class SearchParameters(object):
    def __init__(self):
        self.params = [
            {
                "name": "Stimulus",
                "type": "group",
                "children": [
                    {
                        "name": "Protocol",
                        "type": "list",
                        "values": ["Noise Search", "Tone Search", "Click Search"],
                        "value": "Noise Search",
                    },
                    {
                        "name": "Tone Frequency",
                        "type": "float",
                        "value": 4.0,
                        "step": 1.0,
                        "limits": [0.5, 99.0],
                        "suffix": "kHz",
                        "default": 1.0,
                    },
                    {
                        "name": "Attenuator",
                        "type": "float",
                        "value": 20,
                        "step": 5.0,
                        "limits": [0.0, 120.0],
                        "suffix": "dB",
                        "default": 20.0,
                    },
                    {
                        "name": "Rise-Fall",
                        "type": "float",
                        "value": 2.5,
                        "step": 0.5,
                        "limits": [0.5, 20.0],
                        "suffix": "ms",
                        "default": 2.5,
                    },
                    {
                        "name": "Interstimulus Interval",
                        "type": "float",
                        "value": 1.0,
                        "step": 0.1,
                        "limits": [0.1, 30.0],
                        "suffix": "s",
                        "default": 1.0,
                        "tip": "Time between stimuli in a sweep",
                    },
                    {
                        "name": "Duration",
                        "type": "float",
                        "value": 0.1,
                        "step": 0.05,
                        "limits": [0.01, 10.0],
                        "suffix": "s",
                        "default": 1.0,
                        "tip": "Sound duration, in seconds",
                    },
                    {

                        "name": "Click Duration",
                        "type": "float",
                        "value":   100.,
                        "step":    20.,
                        "limits": [10.0, 1000.0],
                        "suffix": "us",
                        "default": 100.0,
                        "tip": "Click duration, in microseconds",
                    },
                    {
                        "name": "Delay",
                        "type": "float",
                        "value": 0.1,
                        "step": 0.01,
                        "limits": [0.0, 10.0],
                        "suffix": "s",
                        "default": 0.1,
                        "tip": "Delay to stimulus onset from trigger, in seconds",
                    },
                    {
                        "name": "Max Repetitions",
                        "type": "int",
                        "value": 10000,
                        "step": 5,
                        "limits": [1, 10000],
                        "default": 100,
                        "tip": "Max # repetitions to search timeout",
                    },

                ],
            },
        ]

###############################################################################
# Build the window and the GUI controls
###############################################################################

class BuildGui():
    def __init__(self):
        self.app = pg.mkQApp()
        self.mainwin = QtGui.QMainWindow()
        self.win = QtGui.QWidget()
        self.layout = QtGui.QGridLayout()
        self.win.setLayout(self.layout)
        self.mainwin.setCentralWidget(self.win)
        self.mainwin.show()
        self.mainwin.setWindowTitle("Stim Controller")
        self.mainwin.setGeometry(100, 100, 1024, 800)
        self.statusBar = QtGui.QStatusBar()
        self.mainwin.setStatusBar(self.statusBar)
        self.statusMessage = QtGui.QLabel("")
        self.statusBar.addWidget(self.statusMessage)
        self.permStatusMessage = QtGui.QLabel('<b><font color="#00FF00">Ready</b>')
        self.statusBar.addPermanentWidget(self.permStatusMessage)
        self.thread_manager = QThreadPool()  # threading control for using pystim

        self.img = None
        self.timer = QTimer()

        # Define parameters that control aquisition for different paradigms
        SP = SearchParameters()  # the search parameters

        self.ptree = ParameterTree()
        self.ptreedata = Parameter.create(name="params", type="group", children=SP.params)
        self.ptree.setParameters(self.ptreedata)
        ptreewidth = 120

        # now build the ui
        # These are the hardwired buttons
        self.btn_waveform = QtGui.QPushButton("Wave")
        self.btn_spectrum = QtGui.QPushButton("Spectrum")
        self.btn_run_preview = QtGui.QPushButton("Search")
        self.btn_run_pause = QtGui.QPushButton("Pause")
        self.btn_run_record = QtGui.QPushButton("Record")
        self.btn_stop = QtGui.QPushButton("Stop")
        self.btn_quit = QtGui.QPushButton("Quit")
        self.btn_reload = QtGui.QPushButton("Reload")
        self.label_status = QtGui.QLabel("Stopped")
        self.label_status.sizeHint = QtCore.QSize(100, 20)
        self.label_status.setStyleSheet("font:bold 14px; color: red")
        self.stimulus_counter = QtGui.QLabel("0")
        self.stimulus_counter.sizeHint = QtCore.QSize(100, 20)
        self.stimulus_counter.setStyleSheet("font:bold 14px; color:blue;")
        
        hbox = QtGui.QGridLayout()
        hbox.setColumnStretch(0, 1)
        hbox.setColumnStretch(1, 1)
        hbox.setColumnStretch(2, 1)
        hbox.setColumnStretch(3, 1)

        hbox.addWidget(self.btn_quit,        0, 0, 1, 1)
        hbox.addWidget(self.btn_reload,      1, 0, 1, 1)
        hbox.addWidget(self.btn_waveform,    0, 1, 1, 1)
        hbox.addWidget(self.btn_spectrum,    1, 1, 1, 1)
        hbox.addWidget(self.btn_run_preview, 0, 2, 1, 1)
        hbox.addWidget(self.btn_run_pause,   1, 2, 1, 1)
        hbox.addWidget(self.btn_run_record,  2, 2, 1, 1)
        hbox.addWidget(self.btn_stop,        3, 2, 1, 1)
        hbox.addWidget(self.label_status,    3, 0, 1, 1)
        hbox.addWidget(self.stimulus_counter, 3, 1, 1, 1)

        self.layout.addLayout(hbox, 0, 0, 1, 2)

        # build layout for plots and parameters
        self.layout.addWidget(self.ptree, 1, 0, 4, 2)  # Parameter Tree on left
        self.layout.setColumnMinimumWidth(0, ptreewidth)

        # add space for the graphs
        view = pg.GraphicsView()
        glayout = pg.GraphicsLayout(border=(50, 50, 50))
        view.setCentralItem(glayout)
        self.layout.addWidget(view, 0, 2, 5, 3)  # data plots on right
        self.plots = {}
        self.plots["Wave"] = glayout.addPlot()
        self.plots["Wave"].getAxis("left").setLabel("V", color="#ff0000")
        self.plots["Wave"].setTitle("Waveform", color="#ff0000")
        self.plots["Wave"].getAxis("bottom").setLabel("t (sec)", color="#ff0000")
        self.plots["Wave"].setYRange(-1, 1)

        glayout.nextRow()
        self.plots["LongTermSpec"] = glayout.addPlot()
        self.plots["LongTermSpec"].getAxis("left").setLabel("V", color="#ff0000")
        self.plots["LongTermSpec"].setTitle("LongTerm Spectrum", color="#ff0000")
        self.plots["LongTermSpec"].getAxis("bottom").setLabel("F (Hz)", color="#ff0000")

        # Initialize the controller, set the initial parameters, and
        # connect actions and responses to parameter changes

        self.controller = Controller(
            self.ptreedata, self.plots, self.img, self
        )  # we pass a reference to this GUI also, so we can update

        self.controller.getAllParameters(SP.params)  # synchronize parameters with the tree
        self.ptreedata.sigTreeStateChanged.connect(
            self.controller.change
        )  # connect parameters to their updates

        # now connect the buttons
        # self.recentpath = ''
        # -----TFR-----buttons we need
        self.btn_waveform.clicked.connect(self.controller.show_wave)
        self.btn_spectrum.clicked.connect(self.controller.show_spectrogram)
        self.btn_run_preview.pressed.connect(self.controller.start_run_preview)
        self.btn_run_pause.pressed.connect(self.controller.stop_stim)
        self.btn_run_record.pressed.connect(self.controller.start_run_record)
        self.btn_stop.pressed.connect(self.controller.stop_run)
        self.btn_quit.pressed.connect(self.controller.quit)
        self.btn_reload.pressed.connect(self.reload)
        self.timer.setSingleShot(False)  # Make the timer go on forever...
        self.timer.start(20)  # timer times every 20 msec
        self.timer.timeout.connect(self.update_information) # connect to the update method

    def update_information(self):
        """
        The QTimer event loop
        Anything that needs dynamic updating that is not a user-initiated
        action (and is not related to Synapse directly) should go here. 
        Button actions are handled with their own connect callbacks.
        """
        self.controller.show_stimulus_count()

    def reload(self):
        importlib.reload(sound)
        importlib.reload(Utility)

def main():
    gui = BuildGui()

    ## Start Qt event loop unless running in interactive mode.
    ## Event loop will wait for the GUI to activate the updater and start sampling.
    if (sys.flags.interactive != 1) or not hasattr(QtCore, "PYQT_VERSION"):
        QtGui.QApplication.instance().exec_()
    gui.controller.quit()


if __name__ == "__main__":
    main()
