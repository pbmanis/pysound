"""
Stimulus controller
Generates waveforms, controls stimulus presentation
Reiles on sound.py for stimulus waveform generation
Relies on pystim.py for hardware interactions.

Operates in two modes for output intensity:

Atten mode: the levels refer to attenuation in Db, no correction for stimlus
SPL mode : the levels are corrected by the system calibration on a per frequency basis.
    Cannot currently be done for broadband stimuli.

Currently, only the atten mode is supported. 

pbm 2021-2024, tfr 2019-2020
"""

import datetime
import os
from pathlib import Path
import pickle
import pprint
import sys
import time
from collections import OrderedDict

import numpy as np

# from backports import configparser
import pyqtgraph as pg
import pyqtgraph.dockarea
import scipy.io.wavfile as wav
import scipy.signal

# import TDTTankInterface as TDT
import tdt  # necessary for interacting with the tanks in Synapse
from pyqtgraph import QtCore, QtGui, QtWidgets
from pyqtgraph.parametertree import Parameter, ParameterTree

from pysound import Utility, pystim, sound

pp = pprint.PrettyPrinter(indent=4)


class Controller(object):
    def __init__(self, ptreedata, plots, img, maingui):
        self.PS = pystim.PyStim(
            required_hardware=["PA5", "NIDAQ", "RZ5D"],
            controller=self,
        )  # , 'RP21'])
        self.ptreedata = ptreedata
        self.plots = plots  # access to plotting area
        self.img = img
        self.searchmode = False
        # set up a timer to control timing of stimuli
        self.TrialTimer = QtCore.QTimer()  # get a Q timer
        self.TrialTimer.timeout.connect(self.next_stimulus)
        self.maingui = maingui
        self.ProtocolNumber = 0
        self.setAllParameters(ptreedata)
        self.wave = np.zeros(10)

        # Special for clickable map
        # we don't save the data so we don't directly program these -
        # they change with the point clicked in the map

        self.attn = 35
        self.tone_frequency = 4.0  # khz

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
                elif par.type() == 'bool':
                    self.CPars[ch.name()][par.name()] = par.value()

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
            # update grid? 
            if path[0] == 'Stimulus':
                if path[1] == 'Intensities':
                # self.CPars['Stimulus']['Intensities'] = data
                    self.show_FRA_grid(Utility.seqparse(data)[0][0], Utility.seqparse(self.CPars['Stimulus']['Frequencies'])[0][0], clear=True)
                if path[1] == "Frequencies":
                    # self.CPars['Stimulus']['Frequencies'] = data
                    self.show_FRA_grid(Utility.seqparse(self.CPars['Stimulus']['Intensities'])[0][0], Utility.seqparse(data)[0][0], clear=True)

        # self.showParameters()

    def showParameters(self):
        for k in list(self.CPars.keys()):
            print(("Group: %s" % k))
            if isinstance(self.CPars[k], dict):
                for d in list(self.CPars[k].keys()):
                    print(("   %s = %s" % (d, str(self.CPars[k][d]))))
                # pp.pprint(self.CPars)

    # def getCurrentBlock(self):
        # if self.maingui.TT.available:
        #     self.maingui.TT.open_tank()
        #     lastblock = self.maingui.TT.find_last_block()
        #     print ('Current block: ', lastblock)
        #     self.StimRecord['FirstBlock'] = lastblock
        #     self.maingui.TT.close_tank()
        # else:
        #     self.StimRecord['FirstBlock'] = 1

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
        # self.showParameters()
        t = np.linspace(0, 0.1, 1000)
        self.maingui.online_plot.setData([0,0], [0,0]) # , pen=pg.mkPen("g", width=0.30))

        # check for valid times:
        repetitions = self.CPars["Stimulus"]["Repetitions"]  # repetitions in a single sweep
        isi = self.CPars["Stimulus"][
            "InterStimulus Interval"
        ]  # time between stimuli, when doing repeats (repetitions > 1)
       
        self.StimRecord = {}

        # if isi > iti:
        #     msgbox = pg.QtWidgets.QMessageBox()
        #     msgbox.setSizePolicy(
        #         pg.QtWidgets.QSizePolicy.Policy.Expanding,pg.QtWidgets.QSizePolicy.Policy.Expanding
        #     )
        #     msg = "<b><fontcolor: 0xff0000> Stimulus: InterTrial Interval must be >= InterStimulus Interval</b>"
        #     msgbox.setText(msg)
        #     msgbox.exec()
        #     self.maingui.permStatusMessage.setText(msg)
        #     return

        self.runtime = 0
        self.NSamples = 0
        self.running = True
        self.stop_hit = False
        self.startTime = datetime.datetime.now()
        self.StimRecord["StartTime"] = self.startTime  # capture time
        self.StimRecord["Params"] = self.CPars  # all selected parameters
        self.StimRecord["Trials"] = []  # store trial info in a list
        self.StimRecord["savedata"] = (
            True  # flag to save data - set to false by search modes
        )
        self.streamer = None

        # if self.CPars["Stimulus"]["Randomize"]:

        # else:

        # self.getCurrentBlock()
        self.prepare_run()  # reset the data arrays and calculate the next stimulus
        self.lastfreq = None
        self.lastspl = None
        self.trial_count = 0
        self.psth_data = []
        self.isih_data = []
        self.RI_data = []  # dict of attn: [spike counts] (by trials)
        self.RI_levels = []

        self.TrialTimer.setSingleShot(True)
        self.trial_active = True
        self.maingui.label_status.setText("Running")
        self.maingui.label_trialctr.setText("Trial: %04d" % 0)
        self.PS.enable_digital_display=False
        # if "oScope1" in self.PS.RZ5DParams["device_names"]:
        #     self.streamer = tdt.APIStreamer(gizmo='APIStreamer1Ch1', history_seconds=5, callback=self.show_stream_data)
        # print("RZ5D Devices: ", self.PS.RZ5DParams['device_names'])
        gizmonames =self.PS.RZ5D.getGizmoNames()
        # print("Gizmos with API parameters: ",gizmonames)
        streamer_name = "APIStreamerMC1"
        if streamer_name in gizmonames:
            streamer_pars = self.PS.RZ5D.getParameterNames(streamer_name)
            # print("streamer paremters: ", streamer_pars)
            self.PS.RZ5D.getGizmoInfo(streamer_name)
        print("Gizmo Names: ", gizmonames)
        if "APIStreamerMC1" in self.PS.RZ5DParams["device_names"]:
            self.streamer = tdt.APIStreamer(gizmo = streamer_name, history_seconds=2,
                                             callback=self.show_stream_data, verbose=True)
        # print("self.streamer functions: ", dir(self.streamer))
        # print("self.streamer.callback: ", self.streamer.callback)
        else:
            self.PS.enable_digital_display=True
            self.PS.analysis_plots = self.maingui.plots
            self.PS.psth_data = self.psth_data
            self.PS.isih_data = self.isih_data
            self.PS.RI_data = self.RI_data
            self.PS.RI_levels = self.RI_levels
            self.maingui.plots["PSTH"].clear()
            self.maingui.plots["ISIH"].clear()
            self.maingui.plots["RI_plot"].clear()
            self.psth_data = []
            self.isih_data = []

        self.TrialTimer.start(10)  # start (almost) right away  - time is in msec

    def pause_run(self):
        """
        Pause the run - can continue later. This just stops the timer.
        Data is not written until stopping conditions are encountered, or stop is clicked.
        """
        self.maingui.label_status.setText("Paused")
        self.pause_hit = True
        self.TrialTimer.stop()

    def continue_run(self):
        """
        Continue the run if it has been paused.
        """
        if self.pause_hit:
            self.maingui.label_status.setText("Running")
            if self.trial_active:
                self.TrialTimer.start(
                    0.1
                )  # start (almost) right away where we left off
        else:
            return

    def show_stream_data(self, result):
        #  CALLBACK
        # get and display the most recent data
        # print("Streamer Callback show_stream_data called")
        # if self.streamer is None:
        #     print("Streamer is None")
        #     return 
        # print("result: ", dir(result))  # structure
        # print("result.result: ", dir(result.result))  # list
        # print("data lock: ", result.data_lock)
        # print("data shape: ", result.data.shape)
        # print("result data: ", result.data)
        # print("result time array: ", result.time_array)  # not time - but sequential... 
        # print("new data: ", result.new_data.shape)
        ts = result.ts
        # print("ts: ", ts.shape, np.min(ts), np.max(ts))
        if result.data.shape[0] == 0:
            return
        else:
            # print("A..", end="")
            pc = ['g', 'y', 'c', 'b']
            # print("result.data shape 0: ", result.data.shape[0])
            # self.maingui.online_plot.clear()
            for ich in range(result.data.shape[0]):
                if ich != 2:
                    continue
                # print("ich: ", ich, end="")
                # print(result.data[ich,:20])
                # print(result.ts.shape, result.data.shape)
                self.maingui.online_plot.setData(result.ts, result.data[ich,:], pen=pg.mkPen(pc[ich], width=0.30))
                # self.maingui.online_plot.setXRange([np.min(result.ts), np.max(result.ts)])
                 # print("B")

    def next_stimulus(self):
        """
        Present the next stimulus in the sequence

        Parameters
        ----------
        None

        Returns
        -------
        Nothing

        """

      
        self.TrialTimer.stop()
        if self.trial_count > self.total_trials:
            print("Stopped by trial_counter")
            self.stop_run()
            return
        self.maingui.label_trialctr.setText(
            f"Trial: {(self.trial_count + 1):04d} of {self.total_trials:04d}"
        )

        self.TrialTimer.start(
            int(1000.0 * self.CPars["Stimulus"]["InterStimulus Interval"])
        )  # reinit timer
        # do diferently according to protocol:
        spl = self.CPars["Stimulus"]["Attenuator"]
        freq = self.CPars["Stimulus"]["Tone Frequency"]
        protocol = self.CPars["Stimulus"]["Protocol"]
        self.StimRecord["Trials"].append(
            {"Time": "{:%Y.%m.%d %H:%M:%S}".format(datetime.datetime.now())}
        )  # save start time for each trial
        # if self.maingui.TT.available:
        #     self.maingui.TT.open_tank()
        #     lastblock = self.maingui.TT.find_last_block()
        #     self.maingui.TT.close_tank()
        #     self.StimRecord['Trials'][-1]['Block'] = lastblock
        # else:
        #     self.StimRecord['Trials'][-1]['Block'] = 1
        match protocol:
            case "Noise Search" | "Tone Search" |  "Click Search":
                self.play_stimulus(spl=spl, protocol=protocol, frequency=freq, save=False)

            case "One Tone": 
                self.play_stimulus(spl=spl, protocol=protocol,  frequency=freq, save=False)

            case "Tone RI" |"Noise RI":
                spl = self.stim_vary["Intensity"][self.trial_count]
                freq = self.CPars["Stimulus"]["Tone Frequency"]
                print('spl:', spl, "trial: ", self.trial_count)
                print('Protocol {0:s}  attn: {1:3.1f}'.format(protocol, spl))
                self.play_stimulus(protocol=protocol, spl=spl, frequency=freq, save=True)

            case "FRA":
                spl = self.stim_vary["Intensity"][self.trial_count]
                if self.lastspl is None or spl != self.lastspl:
                    time.sleep(self.CPars["Stimulus"]["InterStimulus Interval"])
                freq = self.stim_vary["Frequency"][self.trial_count]
                self.maingui.plots["FRA"].plot([freq, freq],  # change square color
                                                   [spl, spl], symbol='s',
                                                   symbolSize=4, symbolBrush=pg.mkBrush('c'),
                                                   symbolPen=pg.mkPen('c'))
                if (
                    self.lastfreq is None or freq != self.lastfreq
                ):  # determine if we need to calculate the waveform
                    if self.lastfreq is not None:  # add intertrial interval
                        time.sleep(self.CPars["Stimulus"]["InterStimulus Interval"])
                    self.lastfreq = freq
                
                    wave = sound.TonePip(
                        rate=self.PS.Stimulus.out_sampleFreq,
                        duration=self.CPars["Stimulus"]["Duration"]
                        + self.CPars["Stimulus"]["Delay"],
                        f0=freq * 1000.0,
                        dbspl=spl,
                        pip_duration=self.CPars["Stimulus"]["Duration"],
                        pip_starts=[self.CPars["Stimulus"]["Delay"]],
                        ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000,
                    )
                    self.wave = self.map_voltage(protocol, wave.sound, clip=True)
                print(
                    (
                        "Protocol {0:s}  freq: {1:6.3f}  spl: {2:3.1f}".format(
                            protocol, freq, spl
                        )
                    )
                )
                self.play_stimulus(protocol=protocol, spl=spl, frequency=freq, save=True)

            case _:  # all other stimulus sets
                spl = self.CPars["Stimulus"]["Attenuator"]
                self.play_stimulus(protocol=protocol, spl=spl, frequency=None, save=True)

        self.StimRecord["Trials"][-1]["protocol"] = protocol
        self.StimRecord["Trials"][-1]["spl"] = spl
        self.StimRecord["Trials"][-1]["freq"] = freq
        self.trial_count = self.trial_count + 1
        if self.trial_count >= self.total_trials:
            print("Stopping by trial_count >= self.total_trials")
            self.stop_run()  # end last trial without waiting for the rest
            self.PS.stop_recording()
        time.sleep(0.2)  # allow other events

    def play_stimulus(self, protocol:str, spl: float, frequency: float, save:bool=False):
        """
        present stimuli stimuli, at a designated SPL
        """
        self.StimRecord["savedata"] = save
        self.PS.play_sound(
            self.wave,
            self.wave,
            samplefreq=self.PS.Stimulus.out_sampleFreq,
            interstimulus_interval=self.CPars["Stimulus"]["InterStimulus Interval"],
            repetitions=self.CPars["Stimulus"]["Repetitions"],
            attns=self.convert_spl_attn(spl),
            freq=frequency,
            protocol=protocol,
            storedata=self.StimRecord["savedata"],
        )

    def stop_run(self):
        """
        End a run, stop recording, and write the data

        Parameters
        ----------
        None

        Returns
        -------
        Nothing
        """
        self.PS.stop_recording()
        # if self.streamer is not None:
        #     self.streamer.reset()
        self.TrialTimer.stop()

        self.maingui.label_status.setText("Stopped")
        self.trial_active = False
        self.stop_hit = True
        if self.stop_hit == True and self.searchmode == False:
            self.storeData()
        else:  # reset the attenuator and get ready for the next stimulus.
            self.PS.HwOff()
            return

    def storeData(self):
        """
        Write the stimulus parameters to a disk file (ultimately in the current tank)
        """
        # print('first block storedata: ', self.StimRecord['FirstBlock'])
        alldat = [self.CPars, self.StimRecord]
        print(("searchmode:", self.searchmode))

        if (
            self.searchmode is False
        ):  # TFR 20180227- only write the .p file if we're not in search mode
            # desired sequence of events: determine the tank, write the stimulus info into the Tank/Block?
            # TankLocus= self.tdt.SynapseAPI.getCurrentTank()# if self.maingui.TT.available:
            # # TankLocus.replace("\\",'\')
            # if TankLocus!='':
            # print('protocol?: ',self.StimRecord['Trials'][0]['protocol'])
            # return
            # print('Stimulus info: ',self.BlockString)
            print("stim record trials: ", self.StimRecord["Trials"])

            self.BlockString = self.StimRecord["Trials"][0]["protocol"]
            self.BlockString.replace("  ", "")
            print("Blockstring: ", self.BlockString)
            print("Tank name: ", self.PS.TankName)
            filename = Path(self.PS.TankName, f"{self.BlockString:s}.p")
            fh = open(filename, "wb")
            # fh = open(os.path.join(self.maingui.TT.tank_directory,
            #         'Protocol_%s_Blocks_%d-%d.p' % (self.CPars['Stimulus']['Protocol'],
            # self.StimRecord['FirstBlock'], self.StimRecord['Trials'][-1]['Block'])), 'w')
            pickle.dump(alldat, fh)
            fh.close()

    def quit(self):
        self.TrialTimer.stop()
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
        knownprotocols = [
            "Noise Search",
            "Tone Search",
            "Click Search",
            "Tone RI",
            "Single Tone",
            "Noise RI",
            "FRA",
            "RSS",
            "DMR",
            "SSN",
            "Tone SAM",
            "Noise SAM",
            "Clicks",
            "FM Sweep",
            "NotchNoise",
            "Noise Bands",
            "One Tone",
            "CMMR",
            "Noise Train",
        ]
        if protocol not in knownprotocols:
            raise ValueError(
                "Protocol not in list we can map for scaling the voltage in map_voltage"
            )
        if protocol.find("Tone") >= 0 or protocol in ["FRA", "RSS", "One Tone"]:
            print("Tone")
            A = self.CPars["Voltage_Scales"]["Tone_V"]
        if protocol in ["Clicks", "Click Search"]:
            A = self.CPars["Voltage_Scales"]["Click_V"]
        if protocol.find("Noise") >= 0:
            print("Noise")
            A = self.CPars["Voltage_Scales"]["Noise_V"]
        if protocol in ["DMR", "SSN", "FM Sweep"]:
            print("other")
            A = 1.0
        if protocol in ["CMMR"]:
            A = 5.0
        if protocol in ["Noise SAM", "Tone SAM"]:
            A = A / 2.0
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
        Fs = self.PS.Stimulus.out_sampleFreq  # sample frequency
        stim = self.CPars["Stimulus"]["Protocol"]
        #        level = None  # level is dbspl normally for models, but set to None for TDT (5V tone reference)
        seed = 32767
        # print('stim: ', stim)
        wave = None
        self.stim_vary = None
        self.total_trials = 1000
        # print(( 'stim: ', stim))
        if stim in ["Clicks", "Click Search"]:
            wave = sound.ClickTrain(
                rate=Fs,
                duration=self.CPars["Stimulus"]["Duration"],
                dbspl=level,
                click_duration=self.CPars["Clicks"]["Duration"],
                click_starts=1e-3
                * np.arange(
                    self.CPars["Stimulus"]["Delay"] * 1000.0,
                    self.CPars["Clicks"]["Interval"]
                    * (self.CPars["Clicks"]["Number"] + 1)
                    + self.CPars["Stimulus"]["Delay"],
                    self.CPars["Clicks"]["Interval"],
                ),
            )

        elif stim in ["Tone RI", "Tone Search", "Single Tone"]:
            if freq is None:
                freq = self.CPars["Stimulus"]["Tone Frequency"] * 1000.0
            if stim in ["Tone RI"]:
                # print("stimController: stim in Tones: ", (list(self.CPars["Stimulus"].keys())))
                self.stim_vary = {
                    "Intensity": Utility.seqparse(
                        self.CPars["Stimulus"]["Intensities"]
                    )[0][0]
                }
                if self.CPars["Stimulus"]["Randomize"]:
                    rng = np.random.default_rng()
                    self.stim_vary["Intensity"] = rng.choice(self.stim_vary["Intensity"], replace=False)
                self.total_trials = len(self.stim_vary["Intensity"])
            wave = sound.TonePip(
                rate=Fs,
                duration=self.CPars["Stimulus"]["Duration"]
                + self.CPars["Stimulus"]["Delay"],
                f0=freq,
                dbspl=level,
                pip_duration=self.CPars["Stimulus"]["Duration"],
                pip_starts=[self.CPars["Stimulus"]["Delay"]],
                ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000,
            )

        elif stim in ["One Tone"]:
            self.total_trials = 1
            wave = sound.TonePip(
                rate=Fs,
                duration=self.CPars["Stimulus"]["Duration"]
                + self.CPars["Stimulus"]["Delay"],
                f0=self.tone_frequency * 1000,
                dbspl=level,
                pip_duration=self.CPars["Stimulus"]["Duration"],
                pip_starts=[self.CPars["Stimulus"]["Delay"]],
                ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000,
            )

        elif stim in ["Tone SAM"]:
            if freq is None:
                freq = self.CPars["Stimulus"]["Tone Frequency"] * 1000.0
            wave = sound.SAMTone(
                rate=Fs,
                duration=self.CPars["Stimulus"]["Duration"]
                + self.CPars["Stimulus"]["Delay"],
                f0=freq,
                dbspl=level,
                pip_duration=self.CPars["Stimulus"]["Duration"],
                pip_starts=[self.CPars["Stimulus"]["Delay"]],
                ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000.0,
                fmod=self.CPars["SAM"]["Modulation Frequency"],
                dmod=self.CPars["SAM"]["Modulation Depth"],
                seed=seed,
            )

        elif stim in ["FM Sweep"]:
            wave = sound.FMSweep(
                rate=Fs,
                duration=self.CPars["FMSweep"]["Duration"],
                dbspl=level,
                start=self.CPars["Stimulus"]["Delay"],
                ramp=self.CPars["FMSweep"]["Ramp Type"],
                freqs=[
                    self.CPars["FMSweep"]["Freq Start"] * 1000.0,
                    self.CPars["FMSweep"]["Freq End"] * 1000.0,
                ],
            )

        elif stim in ["Noise RI", "Noise Search"]:
            if stim in ['Noise RI']:
                self.stim_vary = {
                    "Intensity": Utility.seqparse(self.CPars["Stimulus"]["Intensities"])[0][
                        0
                    ]
                }
                if self.CPars["Stimulus"]["Randomize"]:
                    rng = np.random.default_rng()
                    # print('intensities: ', self.stim_vary["Intensity"])
                    self.stim_vary["Intensity"] = rng.choice(self.stim_vary["Intensity"],
                                                             size=len(self.stim_vary["Intensity"]),
                                                            replace=False)
                    # print(self.stim_vary["Intensity"])

                self.total_trials = len(self.stim_vary["Intensity"])
            wave = sound.NoisePip(
                rate=Fs,
                duration=self.CPars["Stimulus"]["Duration"]
                + self.CPars["Stimulus"]["Delay"]
                + 0.1,
                dbspl=level,
                pip_duration=self.CPars["Stimulus"]["Duration"],
                pip_starts=[self.CPars["Stimulus"]["Delay"]],
                ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000.0,
                # fmod=self.CPars["SAM"]["Modulation Frequency"],
                # dmod=self.CPars["SAM"]["Modulation Depth"],
                seed=seed,
            )
        elif stim in ["Noise Train"]:
            if self.CPars["Noise Train"]["Number"] > 1:
                howmany = (
                    self.CPars["Noise Train"]["Interval"]
                    + self.CPars["Noise Train"]["Duration"]
                ) * (self.CPars["Noise Train"]["Number"] - 1) + self.CPars["Stimulus"][
                    "Delay"
                ] * 1000
                wave = sound.NoisePip(
                    rate=Fs,
                    duration=self.CPars["Stimulus"]["Duration"]
                    + self.CPars["Stimulus"]["Delay"]
                    + 0,
                    dbspl=level,
                    pip_duration=self.CPars["Noise Train"]["Duration"],
                    pip_starts=1e-3
                    * np.arange(
                        self.CPars["Stimulus"]["Delay"] * 1000.0,
                        howmany,
                        self.CPars["Noise Train"]["Interval"],
                    ),
                    ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000.0,
                    # fmod=self.CPars["CMMR"]["Frequency"],
                    # dmod=0.0,
                    seed=seed,
                )
            else:
                wave = sound.NoisePip(
                    rate=Fs,
                    duration=self.CPars["Stimulus"]["Duration"]
                    + self.CPars["Stimulus"]["Delay"]
                    + 0,
                    dbspl=level,
                    pip_duration=self.CPars["Noise Train"]["Duration"],
                    pip_starts=np.arange(self.CPars["Stimulus"]["Delay"], 1),
                    ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000.0,
                    # fmod=self.CPars["CMMR"]["Frequency"],
                    # dmod=0.0,
                    seed=seed,
                )

        # elif stim in ['Noise RI', 'Noise Search']:
        #     if stim in ['Noise RI']:
        #         self.stim_vary = {'Intensity': Utility.seqparse(self.CPars['Stimulus']['Intensities'])[0][0]}
        #         self.total_trials = len(self.stim_vary['Intensity'])
        #     wave = sound.NoisePip(rate=Fs, duration=self.CPars['Stimulus']['Duration']+self.CPars['Stimulus']['Delay']+0.2,
        #                     f0=self.CPars['Stimulus']['Tone Frequency']*1000., dbspl=level,
        #                     pip_duration=self.CPars['Stimulus']['Duration'], pip_start=[self.CPars['Stimulus']['Delay']],
        #                     ramp_duration=self.CPars['Stimulus']['Rise-Fall']/1000.,
        #                     fmod=self.CPars['CMMR']['Frequency'], dmod=0., seed=seed)
        elif stim in ["Noise SAM"]:
            wave = sound.SAMNoise(
                rate=Fs,
                duration=self.CPars["Stimulus"]["Duration"]
                + self.CPars["Stimulus"]["Delay"],
                f0=self.CPars["Stimulus"]["Tone Frequency"] * 1000.0,
                dbspl=level,
                pip_duration=self.CPars["Stimulus"]["Duration"],
                pip_starts=[self.CPars["Stimulus"]["Delay"]],
                ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000.0,
                fmod=self.CPars["SAM"]["Modulation Frequency"],
                dmod=self.CPars["SAM"]["Modulation Depth"],
                seed=seed,
            )

        elif stim in ["FRA"]:  # frequency response area
            try:
                splseq = Utility.seqparse(self.CPars["Stimulus"]["Intensities"])[0][0]
                freqseq = Utility.seqparse(self.CPars["Stimulus"]["Frequencies"])[0][0]
            except:
                print("SPLS: ", self.CPars["Stimulus"]["Intensities"])
                print(Utility.seqparse(self.CPars["Stimulus"]["Intensities"]))
                print("Freqs: ", self.CPars["Stimulus"]["Frequencies"])
                print(Utility.seqparse(self.CPars["Stimulus"]["Frequencies"]))
                raise ValueError("Unable to parse FRA/Map space")
            self.show_FRA_grid(splseq, freqseq, clear=True)
            mat_spl, mat_freq = np.meshgrid(splseq, freqseq)
            self.stim_vary = {
                "Intensity": mat_spl.ravel(),
                "Frequency": mat_freq.ravel(),
            }

            if self.CPars["Stimulus"]["Randomize"]:
                rng = np.random.default_rng()
                # print('intensities: ', self.stim_vary["Intensity"])
                self.stim_vary["Intensity"] = rng.choice(self.stim_vary["Intensity"],
                                                            size=len(self.stim_vary["Intensity"]),
                                                        replace=False)
                self.stim_vary["Frequency"] = rng.choice(self.stim_vary["Frequency"],
                                                            size=len(self.stim_vary["Frequency"]),
                                                        replace=False)
            print("stim vary: ", self.stim_vary)
            self.total_trials = len(mat_spl.ravel())
            self.last_freq = None
            # note that for this one, the tone is computed at time of use

        elif stim in ["Noise Bands"]:
            self.stim_vary = {
                "Intensity": Utility.seqparse(self.CPars["Stimulus"]["Intensities"])[0][
                    0
                ]
            }
            self.total_trials = len(self.stim_vary["Intensity"])
            wave = sound.NoiseBandPip(
                rate=Fs,
                duration=self.CPars["Stimulus"]["Duration"]
                + self.CPars["Stimulus"]["Delay"],
                f0=self.CPars["Stimulus"]["Tone Frequency"] * 1000.0,
                dbspl=level,
                pip_duration=self.CPars["Stimulus"]["Duration"],
                pip_starts=[self.CPars["Stimulus"]["Delay"]],
                ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000.0,
                seed=seed,
                type=self.CPars["Noise Bands"]["Type"],
                noisebw=self.CPars["Noise Bands"]["Noise BW"] * 1000.0,
                notchbw=self.CPars["Noise Bands"]["Notch BW"] * 1000.0,
                centerfreq=self.CPars["Noise Bands"]["CF"] * 1000.0,
            )

        # elif stim in ["DMR"]:
        #     wave = sound.DynamicRipple(rate=Fs, duration=5.0)

        elif stim in [
            "CMMR"
        ]:  # flanking type is "noise" (modulated), or "MultiTone", or "None".
            # flankingPhase is comodulated or codeviant or random (if type is not None)
            # spacing is band spacing in octaves (for flanking bands)
            #
            wave = sound.ComodulationMasking(rate=Fs,
                                    duration=self.CPars["Stimulus"]["Duration"], 
                                    target_f0 = self.CPars["CMMR"]["Target Frequency"] * 1000.0,
                                    masker_f0 = self.CPars["CMMR"]["Masker Frequency"] * 1000.0,
                                    masker_delay=self.CPars["CMMR"]["Masker Delay"],  # 0.1
                                    masker_duration=self.CPars["CMMR"]["Masker Duration"], # 0.3
                                    target_delay=self.CPars["CMMR"]["Target Delay"], # =0.3, 
                                    target_duration =self.CPars["CMMR"]["Target Duration"], # =0.3,
                                    target_spl= self.CPars["CMMR"]["Target SPL"], # 40,
                                    masker_spl= self.CPars["CMMR"]["Masker SPL"], # 40, 
                                    fmod=self.CPars["CMMR"]["Modulation Frequency"], #10.0, 
                                    dmod=self.CPars["CMMR"]["Modulation Depth"], # 100,
                                    ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000.0, # 0.0025,
                                    flanking_type=self.CPars["CMMR"]["CMMR Flanking Type"],
                                    flanking_spacing=self.CPars["CMMR"]["CMMR Flanking Spacing"], # 0.5, octaves
                                    flanking_phase=self.CPars["CMMR"]["CMMR Flanking Phase"], # 
                                    flanking_bands=self.CPars["CMMR"]["CMMR Flanking Bands"], # 3,
                                    output="Signal", # output could also be other values for testing
                                    )
            
            # Original:
            # wave = sound.ComodulationMasking(
            #     rate=Fs,
            #     duration=self.CPars["Stimulus"]["Duration"]
            #     + self.CPars["Stimulus"]["Delay"],
            #     pip_duration=self.CPars["Stimulus"]["Duration"],
            #     pip_starts=[self.CPars["Stimulus"]["Delay"]],
            #     f0=self.CPars["Stimulus"]["Tone Frequency"] * 1000.0,
            #     ramp_duration=self.CPars["Stimulus"]["Rise-Fall"] / 1000.0,
            #     dbspl=level,
            #     fmod=self.CPars["CMMR"]["Frequency"],
            #     dmod=self.CPars["CMMR"]["Depth"],
            #     flanking_type=self.CPars["CMMR"]["CMMR Flanking Type"],
            #     flanking_spacing=self.CPars["CMMR"]["CMMR Flanking Spacing"],
            #     flanking_phase=self.CPars["CMMR"]["CMMR Flanking Phase"],
            #     flanking_bands=self.CPars["CMMR"]["CMMR Flanking Bands"],
            # )

        elif stim in ["SSN"]:  # speech shaped noise
            # read the file:
            fname = "testsentence.wav"
            (rate, sig) = wav.read(fname)
            duration = float(sig.shape[0] - 1) / rate
            wave = sound.SpeechShapedNoise(
                rate=Fs, duration=duration, waveform=sig[:, 0], samplingrate=rate
            )

        elif stim in ["RSS"]:
            wave = sound.RandomSpectrumShape(
                rate=Fs,
                duration=0.5,
                dbspl=level,
                ramp="linear",
                ramp_duration=1e-2,
                f0=self.CPars["RSS Params"]["CF"] * 1000,
                pip_duration=self.CPars["Stimulus"]["Duration"],
                pip_starts=[self.CPars["Stimulus"]["Delay"]],
                amp_group_size=self.CPars["RSS Params"]["Grouping"],
                amp_sd=self.CPars["RSS Params"]["Level SD"],
                spacing=self.CPars["RSS Params"]["Spacing"],
                octaves=self.CPars["RSS Params"]["Octaves"],
            )

        if stim in ["Noise Search", "Tone Search", "Click Search", "One Tone"]:
            self.searchmode = (
                True  # search mode just runs "forever", until a stop is hit
            )
        else:  # set up the tank to record from and track the presentation no
            self.searchmode = False
            # self.BlockString=stim.replace(" ","")+ProtocolNumber)
            # print('Just made BlockString: ',self.BlockString)
            # self.tdt.SynapseAPI.setCurrentBlock(self.BlockString)
            # self.ProtocolNumber = self.ProtocolNumber + 1

        if wave is not None:
            self.wavesound = wave

            print("prepare_run: len wavesound and time", len(self.wavesound.sound), len(self.wavesound.time))
            print("wave generated")
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
        self.plots["Wave"].plot(self.wavesound.time, self.wave)
        self.plots["Wave"].autoRange(padding=0.05)

    def show_FRA_grid(self, intensities, frequencies, clear: bool=True):
        intens = np.sort(np.unique(intensities))
        freqs = np.sort(np.unique(frequencies))
        minf = np.min(freqs)
        maxf = np.max(freqs)
        minspl = np.min(intens)
        maxspl = np.max(intens)
        if clear:
            self.plots["FRA"].clear()
        for i, db in enumerate(intens):
            self.plots["FRA"].plot([minf, maxf], [db, db], pen=pg.mkPen("darkgrey", width=0.5))
        for i, fr in enumerate(freqs):
            self.plots["FRA"].plot([fr, fr], [minspl, maxspl], pen=pg.mkPen("darkgrey", width=0.5))
        self.plots["FRA"].autoRange()
        self.plots["FRA"].setLogMode(x=True)
        self.show_FRA(intensities, frequencies, clear=False)


    def show_FRA(self, intseries, freqseries, clear: bool = False):
        if isinstance(intseries, str):
            yd = Utility.seqparse(intseries)[0][0]
            xd = Utility.seqparse(freqseries)[0][0]
        else:
            yd=intseries
            xd=freqseries
        xd = np.log10(xd)
        spots = []
        xs=[]
        ys=[]
        data=[]
        for i, y in enumerate(yd):
            for j, x in enumerate(xd):
                spots.append(
                    {
                        "pos": (x, y),
                        "size": 7,
                        "pen": {"color": "darkgrey", "width": 0.5, "alpha": 0.5},
                        "brush": pg.mkBrush("b"),
                        "symbol": "s"
                    }
                )
                xs.append(x)
                ys.append(y)
                data.append((x,y))
        # print("spots: ", spots)
        if clear:
            self.plots["FRA"].clear()
        self.spi = pg.ScatterPlotItem(
            size=7, pen=pg.mkPen("k"), brush=pg.mkBrush("b"), symbol="s"
        )
        self.spi.setData(x=xs, y=ys, hoverable=True, data=data)
        self.plots["FRA"].addItem(self.spi)
        self.plots["FRA"].setLogMode(x=True)

        return self.spi

    def show_spectrogram(self):
        """
        Plot the spectrum in the middle graph

        If spectimage is checked in the main gui, also plot the spectrogram
        in a matplotlib window (must be closed to continue; is blocking)
        """
        self.clearErrMsg()
        self.prepare_run()
        Fs = self.PS.Stimulus.out_sampleFreq
        # show the long term spectrum.
        f, Pxx_spec = scipy.signal.periodogram(
            self.wave, Fs
        )  # , window='flattop', nperseg=8192,
        # noverlap=512, scaling='spectrum')
        self.plots["LongTermSpec"].clear()
        self.plots["LongTermSpec"].plot(f[1:], np.sqrt(Pxx_spec)[1:], pen=pg.mkPen("y"))
        # self.plots['LongTermSpec'].setLogMode(x=True, y=False)

        # print((self.maingui.spectimage))
        if self.maingui.spectimage:  # enable spectrogram plot
            import matplotlib.pyplot as mpl

            ax1 = mpl.subplot(211)
            mpl.plot(self.wavesound.time, self.wave)
            axspec = mpl.subplot(212, sharex=ax1)
            Pxx, freqs, bins, im = mpl.specgram(
                self.wave, NFFT=128, Fs=Fs, noverlap=64, pad_to=256
            )
            # logspec.spectrogram(self.wave, Fs)
            mpl.show()

    #           specfreqs, spectime, Sxx = scipy.signal.spectrogram(self.wavesound.sound*self.Vscale, nperseg=int(0.01*Fs), fs=Fs)
    #           thr = 0. # 1e-8
    #           Sxx[Sxx <= thr] = thr
    # 3 probably better to use matplotlib's spectrogram
    # pos = np.array([0., 1., 0.5, 0.25, 0.75])
    # color = np.array([[0,255,255,255], [255,255,0,255], [0,0,0,255], (0, 0, 255, 255), (255, 0, 0, 255)], dtype=np.ubyte)
    # cmap = pg.ColorMap(pos, color)
    # lut = cmap.getLookupTable(0.0, 1.0, 256)
    # # set colormap
    # # print (dir(self.img))
    # # print (dir(self.img.imageItem))
    # self.img.imageItem.setLookupTable(lut)
    # self.img.setLevels([-40,50])
    # self.img.setImage(Sxx.T)

    def clearErrMsg(self):
        """
        Reset the error notificatoin to the standard ready indicator
        """
        self.maingui.permStatusMessage.setText('<b><font color="#00FF00">Ready</b>')


# Build GUI and window


class BuildGui:
    def __init__(self):
        self.app = pg.mkQApp()
        self.mainwin = pg.GraphicsLayoutWidget()
        self.layout = QtWidgets.QGridLayout()
        self.mainwin.setLayout(self.layout)
        # self.mainwin.setCentralWidget(self.win)
        self.mainwin.show()
        self.mainwin.setWindowTitle("Stim Controller")
        self.mainwin.setGeometry(100, 100, 1024, 800)
        self.spectimage = False
        # self.TDTTankDirectory = ''

        # self.TT = TDT.TDTTankInterface()
        # print('self.TT.available: ', self.TT.available)
        self.statusBar = QtWidgets.QStatusBar()
        self.mainwin.setStatusTip("my status tip")
        self.statusMessage = QtWidgets.QLabel("")
        self.statusBar.addWidget(self.statusMessage)
        self.permStatusMessage = QtWidgets.QLabel('<b><font color="#00FF00">Ready</b>')
        self.statusBar.addPermanentWidget(self.permStatusMessage)

        # retrieve recent path
        # self.configfilename = 'config.ini'
        # if not os.path.isfile(self.configfilename):
        #     # create a configuration file
        #     parser=configparser.SafeConfigParser()
        #     # parser = f.SafeConfigParser()
        #     # initialize parser
        #     parser.add_section('TDTTanks')
        #     parser.set('TDTTanks', 'dir', '')
        #     fh = open(self.configfilename, 'w')
        #     parser.write(fh)
        #     fh.close()
        # else:
        #     parser=configparser.SafeConfigParser()
        #     # parser = ConfigParser.SafeConfigParser()
        #     parser.read('config.ini')
        #     self.TDTTankDirectory = parser.get('TDTTanks', 'dir')
        #     print('tankd dir: ', self.TDTTankDirectory)
        # self.tdt.SynapseAPI.tank_directory = self.tdt.SynapseAPI.getCurrentTank()
        # self.tdt.SynapseAPI.blockstatus = self.tdt.SynapseAPI.setCurrentBlock('test')
        # if len(self.TDTTankDirectory) > 0:
        #     self.TT.open_tank()
        #     lastblock = self.TT.find_last_block()
        #     self.TT.close_tank()
        #     self.TT.show_tank_path()

        # Define parameters that control aquisition and buttons...
        params = [
            {
                "name": "Stimulus",
                "type": "group",
                "children": [
                    {
                        "name": "Protocol",
                        "type": "list",
                        "limits": [
                            "Noise Search",
                            "Tone Search",
                            "Click Search",
                            "Tone RI",
                            "Single Tone",
                            "Noise RI",
                            "FRA",
                            "Clicks",
                            "CMMR",
                            "RSS",
                            # "DMR",
                            "SSN",
                            "Tone SAM",
                            "Noise SAM",
                            "FM Sweep",
                            "Noise Bands",
                            "Noise Train",
                        ],
                        "value": "Noise Search",
                    },
                    {
                        "name": "Tone Frequency",
                        "type": "float",
                        "value": 4.0,
                        "step": 1.0,
                        "limits": [0.5, 99.0],
                        "suffix": "kHz",
                        "default": 4.0,
                    },
                    {
                        "name": "Attenuator",
                        "type": "float",
                        "value": 50,
                        "step": 5.0,
                        "limits": [0.0, 120.0],
                        "suffix": "dB",
                        "default": 50.0,
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
                        "name": "Intensities",
                        "type": "str",
                        "value": "90;20/-10",
                        "suffix": "dBAttn",
                        "default": "90;20/-10",
                    },
                    {
                        "name": "Frequencies",
                        "type": "str",
                        "value": "4;48/8l",
                        "suffix": "kHz",
                        "default": "4;48/8l",
                    },
                    {
                        "name": "Repetitions",
                        "type": "int",
                        "value": 1,
                        "limits": [1, 10000],
                        "default": 1,
                        "tip": "Number of Stimuli per sweep",
                    },
                    # {
                    #     "name": "InterTrial Interval",
                    #     "type": "float",
                    #     "value": 1.0,
                    #     "limits": [0.2, 300.0],
                    #     "suffix": "s",
                    #     "default": 1.0,
                    #     "tip": "Time between sweeps (trials) in FRA and RI protocols",
                    # },
                    {
                        "name": "InterStimulus Interval",
                        "type": "float",
                        "value": 1.0,
                        "limits": [0.2, 30.0],
                        "suffix": "s",
                        "default": 1.0,
                        "tip": "Time between repeated stimuli in a sweep",
                    },
                    {
                        "name": "Randomize",
                        "type": "bool",
                        "value": False,
                        "default": False,
                        "tip": "Randomize presentation order in all dimensions",
                    },
                    {
                        "name": "Duration",
                        "type": "float",
                        "value": 0.2,
                        "step": 0.05,
                        "limits": [0.001, 10],
                        "suffix": "s",
                        "default": 0.2,
                        "tip": "Sound duration, in seconds",
                    },
                    {
                        "name": "Delay",
                        "type": "float",
                        "value": 0.1,
                        "step": 0.05,
                        "limits": [0.001, 10.0],
                        "suffix": "s",
                        "default": 0.1,
                        "tip": "Sound delay, in seconds",
                    },
                ],
            },
            {
                "name": "Clicks",
                "type": "group",
                "expanded": False,
                "children": [
                    {
                        "name": "Interval",
                        "type": "float",
                        "value": 50.0,
                        "step": 5.0,
                        "limits": [1.0, 1000.0],
                        "suffix": "ms",
                        "default": 50.0,
                    },
                    {
                        "name": "Number",
                        "type": "int",
                        "value": 1,
                        "step": 1,
                        "limits": [1, 200.0],
                        "default": 1,
                    },
                    {
                        "name": "Duration",
                        "type": "float",
                        "value": 1e-4,
                        "step": 10e-6,
                        "limits": [10e-6, 1e-3],
                        "suffix": "s",
                        "default": 1e-4,
                    },
                ],
            },
            {
                "name": "SAM",
                "type": "group",
                "expanded": False,
                "children": [
                     {
                        "name": "Modulation Frequency",
                        "type": "float",
                        "value": 10.0,
                        "step": 1.0,
                        "limits": [1.0, 400.0],
                        "suffix": "Hz",
                        "default": 10.0,
                    },
                    {
                        "name": "Modulation Depth",
                        "type": "float",
                        "value": 0.0,
                        "step": 5.0,
                        "limits": [0.0, 200.0],
                        "suffix": "%",
                        "default": 0.0,
                    },
                ]
            },
            {
                "name": "FMSweep",
                "type": "group",
                "expanded": False,
                "children": [
                    {
                        "name": "Duration",
                        "type": "float",
                        "value": 0.5,
                        "step": 0.05,
                        "limits": [5e-3, 10],
                        "suffix": "s",
                        "default": 0.5,
                    },
                    {
                        "name": "Ramp Type",
                        "type": "list",
                        "limits": ["linear", "logarithmic"],
                        "value": "linear",
                    },
                    {
                        "name": "Freq Start",
                        "type": "float",
                        "value": 4,
                        "step": 1,
                        "limits": [1.0, 100.0],
                        "default": 4,
                    },
                    {
                        "name": "Freq End",
                        "type": "float",
                        "value": 48,
                        "step": 1,
                        "limits": [1.0, 100.0],
                        "default": 48,
                    },
                ],
            },
            {
                "name": "CMMR",
                "type": "group",
                "expanded": False,
                "children": [
                    {
                        'name': "Presets",
                        "type": "list",
                        "limits": ["Verhey/Winter", "Target", "OFM+Target", "CM", "CD", "None",],
                        "value": "None",
                        
                    },
                    
                    {
                        "name": "Target Frequency",
                        "type": "float",
                        "value": 4.0,
                        "step": 1.0,
                        "limits": [1.0, 80000.0],
                        "suffix": "kHz",
                        "default": 40.0,
                    },
                     {
                        "name": "Target SPL",
                        "type": "float",
                        "value": 75.0,
                        "step": 5,
                        "limits": [0, 90.],
                        "suffix": "dBSPL",
                        "default": 75,
                    },
                    {
                        "name": "Target Delay",
                        "type": "float",
                        "value": 0.3,
                        "step": 0.01,
                        "limits": [0.1, 5.0],
                        "suffix": "s",
                        "default": 0.3,
                    },
                    {
                        "name": "Target Duration",
                        "type": "float",
                        "value": 0.3,
                        "step": 0.01,
                        "limits": [0.1, 5.0],
                        "suffix": "s",
                        "default": 0.3,
                    },
                     {
                        "name": "Masker Frequency",
                        "type": "float",
                        "value": 4.0,
                        "step": 1.0,
                        "limits": [1.0, 80000.0],
                        "suffix": "kHz",
                        "default": 4.0,
                    },
                     {
                        "name": "Masker SPL",
                        "type": "float",
                        "value": 75.0,
                        "step": 5,
                        "limits": [0, 90.],
                        "suffix": "dBSPL",
                        "default": 75,
                    },
                     {
                        "name": "Masker Delay",
                        "type": "float",
                        "value": 0.1,
                        "step": 0.01,
                        "limits": [0.1, 5.0],
                        "suffix": "s",
                        "default": 0.1,
                    },
                    {
                        "name": "Masker Duration",
                        "type": "float",
                        "value": 0.5,
                        "step": 0.01,
                        "limits": [0.1, 5.0],
                        "suffix": "s",
                        "default": 0.5,
                    },
                     {
                        "name": "Modulation Frequency",
                        "type": "float",
                        "value": 10.0,
                        "step": 1.0,
                        "limits": [1.0, 400.0],
                        "suffix": "Hz",
                        "default": 10.0,
                    },
                    {
                        "name": "Modulation Depth",
                        "type": "float",
                        "value": 100.0,
                        "step": 5.0,
                        "limits": [0.0, 200.0],
                        "suffix": "%",
                        "default": 100.0,
                    },
                    {
                        "name": "CMMR Flanking Type",
                        "type": "list",
                        "limits": ["None", "MultiTone", "NBnoise"],
                        "value": "MultiTone",
                    },
                    {
                        "name": "CMMR Flanking Phase",
                        "type": "list",
                        "limits": ["Comodulated", "Codeviant", "Random"],
                        "value": "Comodulated",
                    },
                    {
                        "name": "CMMR Flanking Bands",
                        "type": "int",
                        "value": 3,
                        "step": 1,
                        "limits": [0, 10],
                        "default": 3,
                    },
                    {
                        "name": "CMMR Flanking Spacing",
                        "type": "float",
                        "value": 0.5,
                        "step": 1 / 8.0,
                        "limits": [1 / 16.0, 2.0],
                        "suffix": "octaves",
                        "default": 0.5,
                    },
                ],
            },

            {
                "name": "RSS Params",
                "type": "group",
                "expanded": False,
                "children": [
                    {
                        "name": "CF",
                        "type": "float",
                        "value": 16.0,
                        "step": 2.0,
                        "limits": [1.0, 64.0],
                        "suffix": "kHz",
                        "default": 16.0,
                    },
                    {
                        "name": "Grouping",
                        "type": "int",
                        "value": 8,
                        "step": 1,
                        "limits": [1, 64],
                        "default": 8,
                    },
                    {
                        "name": "Level SD",
                        "type": "float",
                        "value": 12.0,
                        "step": 2.0,
                        "limits": [0.0, 40.0],
                        "suffix": "dB",
                        "default": 12.0,
                    },
                    {
                        "name": "Spacing",
                        "type": "int",
                        "value": 64,
                        "step": 2,
                        "limits": [1, 128],
                        "suffix": "/octave",
                        "default": 64,
                    },
                    {
                        "name": "Octaves",
                        "type": "float",
                        "value": 3.0,
                        "step": 0.5,
                        "limits": [0.5, 16.0],
                        "default": 3.0,
                    },
                ],
            },
            {
                "name": "Noise Bands",
                "type": "group",
                "expanded": False,
                "children": [
                    {
                        "name": "Type",
                        "type": "list",
                        "limits": ["Bandpass", "BP+Notch"],
                        "value": "Bandpass",
                    },
                    {
                        "name": "Notch BW",
                        "type": "float",
                        "value": 1.0,
                        "step": 1.0,
                        "limits": [0.05, 10.0],
                        "suffix": "kHz",
                        "default": 1.0,
                    },
                    {
                        "name": "Noise BW",
                        "type": "float",
                        "value": 10.0,
                        "step": 1.0,
                        "limits": [1.0, 64],
                        "suffix": "kHz",
                        "default": 10.0,
                    },
                    {
                        "name": "CF",
                        "type": "float",
                        "value": 5.0,
                        "limits": [0.1, 40],
                        "suffix": "kHz",
                        "default": 5.0,
                    },
                ],
            },
            # TFR added this parameter to generate a noise train 20180129
            {
                "name": "Noise Train",
                "type": "group",
                "expanded": False,
                "children": [
                    {
                        "name": "Interval",
                        "type": "float",
                        "value": 100.0,
                        "step": 5.0,
                        "limits": [1.0, 2000.0],
                        "suffix": "ms",
                        "default": 100.0,
                    },
                    {
                        "name": "Number",
                        "type": "int",
                        "value": 2,
                        "step": 1,
                        "limits": [1, 200.0],
                        "default": 2,
                    },
                    {
                        "name": "Duration",
                        "type": "float",
                        "value": 0.05,
                        "step": 10e-2,
                        "limits": [0.001, 1],
                        "suffix": "s",
                        "default": 0.05,
                    },
                ],
            },
            {
                "name": "File From Disk",
                "type": "str",
                "value": "test.wav",
                "default": "test.wav",
            },
        ]

        self.ptree = ParameterTree()
        self.ptreedata = Parameter.create(name="params", type="group", children=params)
        self.ptree.setParameters(self.ptreedata)
        ptreewidth = 120

        # print (dir(self.ptreedata))

        #  print(self.ptreedata.childs)
        #  allpars = OrderedDict()
        #  for ch in self.ptreedata.childs:
        #       allpars[ch.name()] = {}
        #       for par in ch.childs:
        #           print(' name: %s ' % par.name()),
        #           if par.type() == 'int':
        #               print(' %d ' % par.value())
        #               allpars[ch.name()][par.name()] = int(par.value())
        #           elif par.type() == 'float':
        #               allpars[ch.name()][par.name()] = float(par.value())
        #               print(' %f ' % par.value())
        #           elif par.type() == 'list':
        #               print(' %s ' % par.value())
        #               allpars[ch.name()][par.name()] = str(par.value())
        # #          #print( dir(par))
        # #          print('  value: ', par.value(), par.type())
        #  print (allpars)
        #  exit(1)
        # now build the ui
        # hardwired buttons
        self.btn_waveform = QtWidgets.QPushButton("Wave")
        self.btn_spectrum = QtWidgets.QPushButton("Spectrum")
        self.btn_run = QtWidgets.QPushButton("Run")
        self.btn_pause = QtWidgets.QPushButton("Pause")
        self.btn_continue = QtWidgets.QPushButton("Continue")
        self.btn_stop = QtWidgets.QPushButton("Stop")
        self.btn_quit = QtWidgets.QPushButton("Quit")
        self.btn_tdt = QtWidgets.QPushButton("TDT Tank")
        self.label_status = QtWidgets.QLabel("Stopped")
        self.label_trialctr = QtWidgets.QLabel("Trial: 0")
        self.label_status.sizeHint = QtCore.QSize(100, 20)
        self.label_trialctr.setAutoFillBackground(True)
        self.label_trialctr.sizeHint = QtCore.QSize(100, 20)
        self.spect_check = QtWidgets.QCheckBox("Spectrogram")
        self.spect_check.setChecked(False)  # just be sure.
        hbox = QtWidgets.QGridLayout()
        hbox.setColumnStretch(0, 1)
        hbox.setColumnStretch(1, 1)
        hbox.setColumnStretch(2, 1)

        hbox.addWidget(self.btn_quit, 0, 0, 1, 1)
        hbox.addWidget(self.label_status, 0, 1)
        hbox.addWidget(self.btn_run, 0, 2, 1, 1)
        hbox.addWidget(self.btn_pause, 0, 3, 1, 1)

        hbox.addWidget(self.label_trialctr, 1, 0, 1, 2)
        hbox.addWidget(self.btn_stop, 1, 2, 1, 1)
        hbox.addWidget(self.btn_continue, 1, 3, 1, 1)

        hbox.addWidget(self.btn_waveform, 2, 0, 1, 1)
        hbox.addWidget(self.btn_spectrum, 2, 1, 1, 1)
        hbox.addWidget(self.btn_tdt, 2, 2, 1, 1)
        hbox.addWidget(self.spect_check, 2, 3, 1, 1)

        self.layout.addLayout(hbox, 0, 0, 1, 2)

        # build layout for plots and parameters
        self.layout.addWidget(self.ptree, 1, 0, 4, 2)  # Parameter Tree on left
        self.layout.setColumnMinimumWidth(0, ptreewidth)

        # add space for the graphs
        view = pg.GraphicsView()

        self.DockArea = pg.dockarea.DockArea()

        self.Dock_Stimuli = pg.dockarea.Dock("Stimulus Plots")
        self.Dock_Responses = pg.dockarea.Dock("Response Plots")
        self.DockArea.addDock(self.Dock_Stimuli)
        self.DockArea.addDock(self.Dock_Responses, 'below', self.Dock_Stimuli)
        self.layout.addWidget(self.DockArea, 0, 2, 5, 3)
        self.plots = {}
        
        # Stimulus waveform/spectrum plots:
        self.plots["LongTermSpec"] = pg.PlotWidget(title="Long Term Spectrum")
        self.plots["LongTermSpec"].getAxis("left").setLabel("V", color="#ff0000")
        self.plots["LongTermSpec"].setTitle("LongTerm Spectrum", color="#ff0000")
        self.plots["LongTermSpec"].getAxis("bottom").setLabel("F (Hz)", color="#ff0000")
        self.Dock_Stimuli.addWidget(self.plots["LongTermSpec"], 1, 0, 4, 1)

        self.plots["Wave"] = pg.PlotWidget(title="Waveform")
        self.plots["Wave"].getAxis("left").setLabel("V", color="#ff0000")
        self.plots["Wave"].setTitle("Waveform", color="#ff0000")
        self.plots["Wave"].getAxis("bottom").setLabel("t (sec)", color="#ff0000")
        self.plots["Wave"].setYRange(-1, 1)
        self.plots["Wave"].setMaximumHeight(250)
        self.Dock_Stimuli.addWidget(self.plots["Wave"], 0, 0, 1, 1)

        # online analysis plots: Spike raster, PSTH, ISI, RI and FRA
        self.plots["PSTH"] = pg.PlotWidget(title="PSTH") # this is a PlotItem
        self.plots["PSTH"].getAxis("left").setLabel("Counts", color="white")
        self.plots["PSTH"].setTitle("PSTH", color="white")
        self.plots["PSTH"].getAxis("bottom").setLabel("Time", color="white")
        self.plots["PSTH"].setMaximumHeight(200)
        self.PSTH_plot = self.plots["PSTH"].plot([0,0], [0,0], pg.mkPen("g", width=0.35), 
                                                 stepMode='left', fillBrush=pg.mkBrush("g"), fillLevel=0)  # put someting in the PlotItem
        self.Dock_Responses.addWidget(self.plots["PSTH"], row=0, col=0, rowspan=1, colspan=1)

        self.plots["OnLine"] = pg.PlotWidget(title="Spike Raster") # this is a PlotItem
        self.plots["OnLine"].getAxis("left").setLabel("Trial", color="white")
        self.plots["OnLine"].setTitle("Online Analysis", color="white")
        self.plots["OnLine"].getAxis("bottom").setLabel("Time", color="white")
        self.online_plot = self.plots["OnLine"].plot([0,0], [0,0], pg.mkPen("g", width=0.35))  # put someting in the PlotItem
        self.Dock_Responses.addWidget(self.plots["OnLine"], row=1, col=0, rowspan=3, colspan=1)

        self.plots["ISIH"] = pg.PlotWidget(title="ISI") # this is a PlotItem
        self.plots["ISIH"].getAxis("left").setLabel("Counts", color="white")
        self.plots["ISIH"].setTitle("ISIH", color="white")
        self.plots["ISIH"].getAxis("bottom").setLabel("Time", color="white")
        self.plots["ISIH"].setMaximumHeight(200)
        self.ISIH_plot = self.plots["ISIH"].plot([0,0], [0,0], pg.mkPen("g", width=0.35), 
                                                 stepMode='left', fillBrush=pg.mkBrush("b"), fillLevel=0)  # put someting in the PlotItem
        self.Dock_Responses.addWidget(self.plots["ISIH"], row=0, col=1, rowspan=1, colspan=1)

        self.plots["RI_plot"] = pg.PlotWidget(title="Rate-Intensity") # this is a PlotItem
        self.plots["RI_plot"].getAxis("left").setLabel("Spike Counts", color="white")
        self.plots["RI_plot"].setTitle("Rate-Intensity", color="white")
        self.plots["RI_plot"].getAxis("bottom").setLabel("Intensity", color="white")
        self.RI_plot = self.plots["RI_plot"].plot([0,0], [0,0], pg.mkPen("g", width=0.35))  # put someting in the PlotItem
        self.Dock_Responses.addWidget(self.plots["RI_plot"], row=1, col=1, rowspan=3, colspan=1)

        #     if self.spectimage:
        #         self.img = pg.ImageView() # view=self.plots['Spec'])
        #         arr = np.random.random((100, 32))
        #         self.img.setImage(arr)
        #         self.img.ui.roiBtn.hide()
        # #        self.img.ui.menuBtn.hide()
        #         self.img.show()
        #     else:
        self.img = None

        # Frequency Response Area
        self.plots["FRA"] = pg.PlotWidget(title="FRA")
        #        self.l2.addWidget(self.plots["FRA"])
        self.plots["FRA"].getAxis("bottom").setLabel("F (kHz)")
        self.plots["FRA"].getAxis("left").setLabel("dB ATTN")
        self.plots["FRA"].setTitle("FRA")
        self.plots["FRA"].setXRange(0, 50, padding=0)
        # self.plots["FRA"].setLogMode(x=True)
        self.plots["FRA"].setYRange(125, -5, padding=0)
        self.Dock_Responses.addWidget(self.plots["FRA"], row=4, col=0, rowspan=2, colspan=2)

        #     xd = np.arange(2, 48, 1)
        #    # xd = np.logspace(np.log2(2), np.log2(64), 50, base=2)
        #    # print ('xd: ', xd)
        #     yd = np.arange(120, 5, -5)
        #     spots = []
        #     self.lastPoint = None
        #     for i in range(xd.shape[0]):
        #         for j in range(yd.shape[0]):
        #             spots.append({'pos': (xd[i], yd[j]), 'size': 7, 'pen': {'color': 'k', 'width': 0.5, 'alpha': 0.5},
        #                 'brush': pg.mkBrush('b')})
        #     self.spi = pg.ScatterPlotItem(size=7, pen=pg.mkPen('k'), brush=pg.mkBrush('b'), symbol='s')
        # self.spi.addPoints(spots)
        # self.plots["FRA"].addItem(self.spi)
        # self.spi.getViewBox().invertY(True)
        # self.spi.sigClicked.connect(self.getClickedLocation)
        # cross hair
        # vLine = pg.InfiniteLine(angle=90, movable=True)
        # hLine = pg.InfiniteLine(angle=0, movable=True)
        # self.plots["FRA"].addItem(vLine, ignoreBounds=False)
        # self.plots["FRA"].addItem(hLine, ignoreBounds=False)
        # vb = self.plots["FRA"].vb

        # def mouseMoved(evt):
        #     pos = evt[0]  ## using signal proxy turns original arguments into a tuple
        #     if self.plots["FRA"].sceneBoundingRect().contains(pos):
        #         mousePoint = vb.mapSceneToView(pos)
        #         index = int(mousePoint.x())
        #         if index > 0 and index < len(data1):
        #             label.setText("<span style='font-size: 12pt'>x=%0.1f,   <span style='color: red'>y1=%0.1f</span>,   <span style='color: green'>y2=%0.1f</span>" % (mousePoint.x(), data1[index], data2[index]))
        #         vLine.setPos(mousePoint.x())
        #         hLine.setPos(mousePoint.y())
        # proxy = pg.SignalProxy(self.plots["FRA"].scene().sigMouseMoved, rateLimit=60, slot=mouseMoved)

        # self.plots['Plot2'] = pg.plot(Title="Plot2")
        #  self.l2.addWidget(self.plots['Plot2'])
        #  self.plots['Plot2'].getAxis('bottom').setLabel('t (s)')
        #  self.plots['Plot2'].getAxis('left').setLabel('V')
        #  self.plots['Plot2'].setTitle('Plot 2')
        #
        #  self.plots['Plot3'] = pg.plot(Title="Plot3")
        #  self.l2.addWidget(self.plots['Plot3'])
        #  self.plots['Plot3'].setTitle('Plot3')
        #  self.plots['Plot3'].getAxis('bottom').setLabel('t (s)')
        #  self.plots['Plot3'].getAxis('left').setLabel('V')

        #
        # Initialize the controller, set parameters, and connect actions and
        # responses to parameter changes
        #
        self.controller = Controller(
            self.ptreedata, self.plots, self.img, self
        )  # we pass the gui also

        self.controller.setAllParameters(params)  # synchronize parameters with the tree
        self.ptreedata.sigTreeStateChanged.connect(
            self.controller.change
        )  # connect parameters to their updates

        # now connect the buttons
        self.recentpath = ""
        self.btn_waveform.clicked.connect(self.controller.show_wave)
        self.btn_spectrum.clicked.connect(self.controller.show_spectrogram)
        # self.btn_tdt.clicked.connect(print('TDT Tank: ', self.tdt.SynapseAPI.getCurrentTank())) # self.TT.set_tank_path)
        #    self.ButtonEvents = QtCore.QTimer() # get a Q timer
        #    self.btn_stop.clicked.connect(timer, SIGNAL(timeout()), this, SLOT(processOneThing()));
        # timer->start();

        self.btn_run.clicked.connect(self.controller.start_run)

        self.btn_pause.clicked.connect(self.controller.pause_run)
        self.btn_continue.clicked.connect(self.controller.next_stimulus)
        self.btn_stop.clicked.connect(self.controller.stop_run)
        self.btn_quit.clicked.connect(self.controller.quit)
        self.spect_check.clicked.connect(self.speccheck)
        # update the fra plot
        self.controller.show_FRA_grid(Utility.seqparse(self.controller.CPars["Stimulus"]["Intensities"])[0][0],
                            Utility.seqparse(self.controller.CPars["Stimulus"]["Frequencies"])[0][0], 
                             clear=True)
        self.spi = self.controller.show_FRA(
            self.controller.CPars["Stimulus"]["Intensities"],
            self.controller.CPars["Stimulus"]["Frequencies"],
            clear=False,
        )  # first time through, get self.spi.
        self.spi.getViewBox().invertY(True)
        self.spi.sigClicked.connect(self.getClickedLocation)
        # self.updateStatusMessage()

    def speccheck(self):
        self.spectimage = self.spect_check.isChecked()

    # def getTDTTank(self, dirname=None):
    #     filedialog = QtGui.QFileDialog()
    #     filedialog.setFileMode(QtGui.QFileDialog.Directory)

    #         self.TT.tank_directory = str(filedialog.getExistingDirectory(None, "Select Tank Directory", self.recentpath,
    #                                     QtGui.QFileDialog.ShowDirsOnly))
    #         self.recentpath = self.TT.tank_directory
    # #        print('Tank dir selected: ', self.TT.tank_directory)
    #         self.setTankIni(self.TT.tank_directory)
    #         self.TT.open_tank()
    #         lastblock = self.TT.find_last_block()
    #         self.TT.close_tank()
    #         self.TT.show_tank_path()
    #         self.updateStatusMessage()

    # def updateStatusMessage(self):
    #     if self.TT.available is False:
    #         message = ('No TDT Tank')
    #     else:
    #         message = ('Tank: {0:s}  CurrentBlock: {1:d}'.format(self.TT.tank_directory, self.TT.lastblock))
    #     self.statusMessage.setText(message)

    def getClickedLocation(self, points):
        # print (dir(points))
        # print (points.event())
        # print('mouse click: ', points.mouseClickEvent(points.event()))
        # print('mouse doubleclick: ', points.mouseDoubleClickEvent(points.event()))
        self.mousePoint = points.ptsClicked[0].viewPos()
        print(("mousepoint: ", self.mousePoint.x(), self.mousePoint.y()))
        points.ptsClicked[0].setBrush(pg.mkBrush("r"))
        points.ptsClicked[0].setSize(7)
        if self.lastPoint is None:
            self.lastPoint = points.ptsClicked[0]
        else:
            self.lastPoint.setBrush(pg.mkBrush("b"))
            self.lastPoint.setSize(7)
            self.lastPoint = points.ptsClicked[0]

        stimpars = list(self.ptreedata.param("Stimulus").items.keys())[
            0
        ]  # force to One Tone mode
        stimpars.param.names["Protocol"].setValue("One Tone")
        #        stimpars.param.emitStateChanged()  # trigger
        print("One Tone for me!")

        self.controller.protocol = stimpars.param.names["Protocol"].value()
        self.controller.tone_frequency = self.mousePoint.x()
        self.controller.attn = self.mousePoint.y()
        self.searchmode = True
        self.controller.CPars["Stimulus"]["Protocol"] = "One Tone"
        self.controller.prepare_run(
            freq=self.controller.tone_frequency, level=self.controller.attn
        )
        self.controller.start_run()


def main():
    gui = BuildGui()

    ## Start Qt event loop unless running in interactive mode.
    ## Event loop will wait for the GUI to activate the updater and start sampling.
    if (sys.flags.interactive != 1) or not hasattr(QtCore, "PYQT_VERSION"):
        QtGui.QGuiApplication.instance().exec()
    gui.controller.quit()


if __name__ == "__main__":
    main()
