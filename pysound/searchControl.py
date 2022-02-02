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
from dataclasses import dataclass
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


class Controller(object):
    def __init__(self, ptreedata, plots, img, maingui):
        self.ptreedata = ptreedata
        self.plots = plots  # access to plotting area
        self.img = img
        self.running = False
        self.NSamples = 0

        self.maingui = maingui
        # self.ProtocolNumber = 0
        self.TDTinfo = tdt.SynapseAPI()
        self.setAllParameters(ptreedata)

        self.PS = pystim.PyStim(
            required_hardware=["PA5", "NIDAQ", "RZ5D"], controller=self,
        )  # , 'RP21'])

        # Special for clickable map
        # we don't save the data so we don't directly program these - they change with the point clicked in the map

    # --------TFR---- we don't need these
    # self.attn = 35
    # self.tone_frequency = 4.0 # khz

    def reload(self):
        importlib.reload(sound)
        importlib.reload(Utility)

    def setAllParameters(self, params):
        """
        Set all of the local parameters from the parameter tree

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

        # self.showParameters()

    def showParameters(self):
        for k in list(self.CPars.keys()):
            print(("Group: %s" % k))
            if isinstance(self.CPars[k], dict):
                for d in list(self.CPars[k].keys()):
                    print(("   %s = %s" % (d, str(self.CPars[k][d]))))

    def stimulus_count(self):
        self.NSamples += 1

    def start_run(self):
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
        self.NSamples = 0
        self.running = True
        self.stop_hit = False
        self.startTime = datetime.datetime.now()
        self.maingui.label_status.setText("Running")
        self.prepare_run()  # reset the data arrays and calculate the next stimulus
        self.running = True
        self.PS.play_sound(
                wavel=self.wave,
                waver=self.wave,
                samplefreq=self.PS.Stimulus.out_sampleFreq,
                isi=self.CPars["Stimulus"]["Interstimulus Interval"],
                reps=self.CPars["Stimulus"]["Max Repetitions"],
                attns=self.convert_spl_attn(self.CPars["Stimulus"]["Attenuator"]),
                storedata=False,
            )

    def stop_run(self):
        """
        End a run, and write data
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Nothing
        """
        self.running = False
        self.PS.stop_stim()
        self.maingui.label_status.setText("Stopped")
        self.NSamples = 0
        # self.trial_active = False
        # self.stop_hit = True
        # if self.stop_hit == True and self.searchmode == False:
        #     self.storeData()
        # # else: #reset the attenuator and get ready for the next stimulus.
        # self.PS.HwOff()
        return

    def quit(self):
        """
        Provide a clean exit, by stopping all stimuli and disengaging hardware
        """
        # self.TrialTimer.stop()
        self.PS.stop_stim()
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
            raise ValueError("intensity mode must be attenuation or spl")

    def map_voltage(self, protocol, wave, clip=True):
        """
        Provide scaling of the stimulus voltage based on the stimuli. 
        General rule: as high a voltage as possible, but also avoiding
        the possiblity of clipping
        
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

    def prepare_run(self, freq=None, level=None):
        """
        Clear out all arrays for the data collection run
        
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
        self.prepare_run()  # force computation/setup of stimulus
        self.plots["Wave"].clear()
        # self.plots['Wave'].plot(self.wavesound.time, self.wave)
        self.plots["Wave"].plot(self.wavesound.time, self.wave)

    def show_spectrogram(self):
        """
        Plot the spectrum in the middle graph
        
        If spectimage is checked in the main gui, also plot the spectrogram
        in a matplotlib window (must be closed to continue; is blocking)
        """
        self.clearErrMsg()
        self.prepare_run()
        # show the long term spectrum.
        f, Pxx_spec = scipy.signal.periodogram(
            self.wave, self.PS.Stimulus.out_sampleFreq
        )  # , window='flattop', nperseg=8192,
        # noverlap=512, scaling='spectrum')
        self.plots["LongTermSpec"].clear()
        self.plots["LongTermSpec"].plot(f[1:], np.sqrt(Pxx_spec)[1:], pen=pg.mkPen("y"))
        # self.plots['LongTermSpec'].setLogMode(x=True, y=False)

    def clearErrMsg(self):
        """
        Reset the error notificatoin to the standard ready indicator
        """
        self.maingui.permStatusMessage.setText('<b><font color="#00FF00">Ready</b>')


# Build GUI and window


class BuildGui:
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

        # Define parameters that control aquisition and buttons...
        params = [
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

        self.ptree = ParameterTree()
        self.ptreedata = Parameter.create(name="params", type="group", children=params)
        self.ptree.setParameters(self.ptreedata)
        ptreewidth = 120

        # now build the ui
        # hardwired buttons
        self.btn_waveform = QtGui.QPushButton("Wave")
        self.btn_spectrum = QtGui.QPushButton("Spectrum")
        self.btn_run = QtGui.QPushButton("Run")
        self.btn_stop = QtGui.QPushButton("Stop")
        self.btn_quit = QtGui.QPushButton("Quit")
        self.btn_reload = QtGui.QPushButton("Reload")
        self.label_status = QtGui.QLabel("Stopped")
        self.label_status.sizeHint = QtCore.QSize(100, 20)
        # self.label_trialctr.setAutoFillBackground(True)
        # self.label_trialctr.sizeHint = QtCore.QSize(100, 20)

        hbox = QtGui.QGridLayout()
        hbox.setColumnStretch(0, 1)
        hbox.setColumnStretch(1, 1)
        hbox.setColumnStretch(2, 1)

        hbox.addWidget(self.btn_quit, 0, 0, 1, 1)
        hbox.addWidget(self.btn_reload, 1, 0, 1, 1)
        hbox.addWidget(self.btn_run, 0, 2, 1, 1)
        hbox.addWidget(self.btn_stop, 1, 2, 1, 1)
        hbox.addWidget(self.btn_waveform, 2, 0, 1, 1)
        hbox.addWidget(self.btn_spectrum, 2, 1, 1, 1)
        hbox.addWidget(self.label_status, 2, 2, 1, 1)

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

        # Initialize the controller, set parameters, and connect actions and
        # responses to parameter changes
        #
        self.controller = Controller(
            self.ptreedata, self.plots, self.img, self
        )  # we pass a reference to the GUI also

        self.controller.setAllParameters(params)  # synchronize parameters with the tree
        self.ptreedata.sigTreeStateChanged.connect(
            self.controller.change
        )  # connect parameters to their updates

        # now connect the buttons
        # self.recentpath = ''
        # -----TFR-----buttons we need
        self.btn_waveform.clicked.connect(self.controller.show_wave)
        self.btn_spectrum.clicked.connect(self.controller.show_spectrogram)
        self.btn_run.pressed.connect(self.controller.start_run)
        self.btn_stop.pressed.connect(self.controller.stop_run)
        self.btn_quit.pressed.connect(self.controller.quit)
        self.btn_reload.pressed.connect(self.controller.reload)
        self.timer.setSingleShot(False)
        self.timer.start(20)  # timer times every 20 msec
        self.timer.timeout.connect(self.showtime)

    def showtime(self):
        # print('showtime!')
        time.sleep(0.002)
        # if not self.controller.running:
        #     print("Not Running")


def main():
    gui = BuildGui()

    ## Start Qt event loop unless running in interactive mode.
    ## Event loop will wait for the GUI to activate the updater and start sampling.
    if (sys.flags.interactive != 1) or not hasattr(QtCore, "PYQT_VERSION"):
        QtGui.QApplication.instance().exec_()
    gui.controller.quit()


if __name__ == "__main__":
    main()
