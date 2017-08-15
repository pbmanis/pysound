from __future__ import print_function
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

import sys
import os
import datetime
import numpy as np
import scipy.signal
import time
import pickle
import scipy.io.wavfile as wav
from collections import OrderedDict
import ConfigParser
import pyqtgraph as pg
from PyQt4 import QtGui, QtCore
from pyqtgraph.parametertree import Parameter, ParameterTree
import pystim
import Utility
import sound
import pprint
import TDTTankInterface as TDT

pp = pprint.PrettyPrinter(indent=4)


class Controller(object):
    def __init__(self, ptreedata, plots, img, maingui):
        self.PS = pystim.PyStim(hdw=['PA5', 'NIDAQ', 'RZ5D'])
        self.ptreedata = ptreedata
        self.plots = plots  # access to plotting area
        self.img = img
        # set up a timer to control timing of stimuli
        self.TrialTimer=QtCore.QTimer() # get a Q timer
        self.TrialTimer.timeout.connect(self.next_stimulus);
        self.maingui = maingui
        
        self.setAllParameters(ptreedata)

        # Special for clickable map
        # we don't save the data so we don't directly program these - they change with the point clicked in the map

        self.attn = 35
        self.tone_frequency = 4.0 # khz 
                
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
        self.CPars['IntensityMode'] = 'attenuation'
        self.CPars['Voltage_Scales'] = {'Tone_V': 10.0, 'maxTone_dB': {'MF1': 110, 'EC1': 83.9},
                   'Click_V': 5.0, 'maxClick_dB': {'MF1': 108.5, 'EC1': 79.5},
                   'Noise_V': 2.5, 'maxNoise_dB': {'MF1': 0, 'EC1': 0}}  # we don't actualy know... but also need to clip
        
        for ch in self.ptreedata.childs:
            self.CPars[ch.name()] = {}
            for par in ch.childs:
                #print(' name: %s ' % par.name()),
                if par.type() == 'int':
                    self.CPars[ch.name()][par.name()] = int(par.value())
                elif par.type() == 'float':
                    self.CPars[ch.name()][par.name()] = float(par.value())
                elif par.type() == 'list':
                    self.CPars[ch.name()][par.name()] = str(par.value())
                elif par.type() == 'str':
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
        
        #self.showParameters()

    def showParameters(self):
        for k in self.CPars.keys():
            print('Group: %s' % k)
            if isinstance(self.CPars[k], dict):
                for d in self.CPars[k].keys():
                    print('   %s = %s' % (d, str(self.CPars[k][d])))
                #pp.pprint(self.CPars)

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
        
        # check for valid times: 
        nr = self.CPars['Stimulus']['Repetitions']  # repetitions in a single sweep
        isi = self.CPars['Stimulus']['Interstimulus Interval']  # time between stimuli, s
        iti = self.CPars['Stimulus']['Intertrial Interval']  # time between trials in RI and FRA, s
        self.StimRecord = {}
        
        sweepdur = nr*isi
        if sweepdur > 0.8*iti:
            self.maingui.permStatusMessage('<b><fontcolor: 0xff0000> Stimuli nreps*isi must be < 80\% of iti</b>')
            return
        
        self.runtime = 0
        self.NSamples = 0
        self.running = True
        self.stop_hit = False
        self.startTime = datetime.datetime.now()
        self.StimRecord['StartTime'] = self.startTime  # capture time
        self.StimRecord['Params'] = self.CPars  # all selected parameters
        self.StimRecord['Trials'] = []  # store trial info in a list
        self.StimRecord['savedata'] = True  # flag to save data - set to false by search modes
        if self.maingui.TT.available:
            self.maingui.TT.open_tank()
            lastblock = self.maingui.TT.find_last_block()
            self.maingui.TT.close_tank()
            self.StimRecord['FirstBlock'] = lastblock
        else:
            self.StimRecord['FirstBlock'] = 1
        self.prepare_run()  # reset the data arrays and calculate the next stimulus
        self.lastfreq = None
        self.trial_count = 0
        self.TrialTimer.setSingleShot(True)
        self.trial_active = True
        self.maingui.label_status.setText('Running')
        self.maingui.label_trialctr.setText('Trial: %04d' % 0)
        
        self.TrialTimer.start(0.1) # start (almost) right away

    def pause_run(self):
        """
        Pause the run - can continue later. This just stops the timer.
        Data is not written until stopping conditions are encountered, or stop is clicked.
        """
        self.maingui.label_status.setText('Paused')
        self.pause_hit = True
        self.TrialTimer.stop()

    def continue_run(self):
        """
        Continue the run if it has been paused. 
        """
        if self.pause_hit:
            self.maingui.label_status.setText('Running')
            if self.trial_active:
                self.TrialTimer.start(0.1) # start (almost) right away where we left off
        else:
            return

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
        if self.trial_count >= self.total_trials:
            self.stop_run()
            return
        self.maingui.label_trialctr.setText('Trial: {0:04d} of {1:04d}'.format(self.trial_count+1, self.total_trials))
        
        self.TrialTimer.start(int(1000.0*self.CPars['Stimulus']['Intertrial Interval']))  # reinit timer
         # do diferently according to protocol:
        spl = self.CPars['Stimulus']['Attenuator']
        freq = self.CPars['Stimulus']['Tone Frequency']
        protocol = self.CPars['Stimulus']['Protocol']
        self.StimRecord['Trials'].append({'Time': '{:%Y.%m.%d %H:%M:%S}'.format(datetime.datetime.now())})  # start time for each trial
        if self.maingui.TT.available:
            self.maingui.TT.open_tank()
            lastblock = self.maingui.TT.find_last_block()
            self.maingui.TT.close_tank()
            self.StimRecord['Trials'][-1]['Block'] = lastblock
        else:
            self.StimRecord['Trials'][-1]['Block'] = 1
        if protocol in ['Noise Search', 'Tone Search']:
            self.StimRecord['savedata'] = False
            self.PS.play_sound(self.wave, self.wave,
                samplefreq=self.PS.out_sampleFreq,
                isi=self.CPars['Stimulus']['Interstimulus Interval'],
                reps=self.CPars['Stimulus']['Repetitions'], 
                attns=self.convert_spl_attn(spl), storedata=self.StimRecord['savedata'])
            
        elif protocol in ['One Tone']:
            spl = self.attn
            self.StimRecord['savedata'] = False
            self.PS.play_sound(self.wave, self.wave,
                samplefreq=self.PS.out_sampleFreq,
                isi=self.CPars['Stimulus']['Interstimulus Interval'],
                reps=self.CPars['Stimulus']['nreps'],
                attns=self.convert_spl_attn(spl), storedata=self.StimRecord['savedata'])
            
        elif protocol in ['Tone RI', 'Noise RI']:
            spl = self.stim_vary['Intensity'][self.trial_count]
            freq = self.CPars['Stimulus']['Tone Frequency']
            # print('spl:', spl)
            # print('Protocol {0:s}  attn: {1:3.1f}'.format(protocol, spl))
            self.PS.play_sound(self.wave, self.wave,
                samplefreq=self.PS.out_sampleFreq,
                isi=self.CPars['Stimulus']['Interstimulus Interval'],
                reps=self.CPars['Stimulus']['Repetitions'], attns=self.convert_spl_attn(spl))

        elif protocol in ['FRA']:
            spl = self.stim_vary['Intensity'][self.trial_count]
            freq = self.stim_vary['Frequency'][self.trial_count]
            if self.lastfreq is None or freq != self.lastfreq:  # determine if we need to calculate the waveform
                self.lastfreq  = freq
                wave = sound.TonePip(rate=self.PS.out_sampleFreq,
                            duration=self.CPars['Stimulus']['Duration']+self.CPars['Stimulus']['Delay'],
                            f0=freq*1000., dbspl=spl, 
                            pip_duration=self.CPars['Stimulus']['Duration'], pip_start=[self.CPars['Stimulus']['Delay']],
                            ramp_duration=self.CPars['Stimulus']['Rise-Fall']/1000)
                self.wave = self.map_voltage(protocol, wave.sound, clip=True)
            print('Protocol {0:s}  freq: {1:6.3f}  spl: {2:3.1f}'.format(protocol, freq, spl))
            self.PS.play_sound(self.wave, self.wave,
                samplefreq=self.PS.out_sampleFreq,
                isi=self.CPars['Stimulus']['Interstimulus Interval'], reps=self.CPars['Stimulus']['Repetitions'],
                attns=self.convert_spl_attn(spl))

        else:
            spl = self.CPars['Stimulus']['Attenuator']
            self.PS.play_sound(self.wave, self.wave,
                samplefreq=self.PS.out_sampleFreq,
                isi=self.CPars['Stimulus']['Interstimulus Interval'],
                reps=self.CPars['Stimulus']['Repetitions'], attns=self.convert_spl_attn(spl))
        
        self.StimRecord['Trials'][-1]['protocol'] = protocol
        self.StimRecord['Trials'][-1]['spl'] = spl
        self.StimRecord['Trials'][-1]['freq'] = freq
        self.trial_count = self.trial_count + 1
        if self.trial_count >= self.total_trials:
            self.stop_run()  # end last trial without waiting for the rest
        time.sleep(0.05)   # allow other events
        
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
        self.TrialTimer.stop()
        self.maingui.label_status.setText('Stopped')
        self.trial_active = False
        self.stop_hit = True
        self.storeData()
        
    def storeData(self):
        """
        Write the stimulus parameters to a disk file (ultimately in the current tank)
        """
        alldat = [self.CPars, self.StimRecord]
        if self.maingui.TT.available:
            fh = open(os.path.join(self.maingui.TT.tank_directory,
                    'Protocol_%s_Blocks_%d-%d.p' % (self.CPars['Stimulus']['Protocol'], 
            self.StimRecord['FirstBlock'], self.StimRecord['Trials'][-1]['Block'])), 'w')
            pickle.dump(alldat, fh)
            fh.close()
        
    def quit(self):
        self.TrialTimer.stop()
        self.PS.HwOff()
        exit(0)
    
    def convert_spl_attn(self, spl):
        if self.CPars['IntensityMode'] == 'attenuation':
            return spl # use as attenuation directly
        elif self.CPars['IntensityMode'] == 'spl':
            return [100.-spl, 100.-spl]  # just rough, within 15dB, need to clean up
        else:
            raise ValueError('intensity mode must be attenuation or spl')

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
        knownprotocols = ['Noise Search', 'Tone Search',
                        'Tone RI', 'Noise RI', 'FRA',
                        'RSS', 'DMR', 'SSN', 'Tone SAM', 'Noise SAM', 'Clicks', 'FM Sweep',
                        'NotchNoise', 'Noise Bands', 'One Tone', 'CMMR']
        if protocol not in  knownprotocols:
            raise ValueError('Protocol not in list we can map for scaling the voltage in map_voltage')
        if protocol.find('Tone') >=0 or protocol in ['FRA', 'RSS', 'One Tone']:
            print ('Tone')
            A = self.CPars['Voltage_Scales']['Tone_V']
        if protocol.find('Clicks') >= 0:
            A = self.CPars['Voltage_Scales']['Click_V']
        if protocol.find('Noise') >= 0:
            print('Noise')
            A = self.CPars['Voltage_Scales']['Noise_V']
        if protocol in ['DMR', 'SSN', 'FM Sweep']:
            print('other')
            A = 1.0
        if protocol in ['CMMR']:
            A = 5.0
        if protocol in ['Noise SAM', 'Tone SAM']:
            A = A / 2.0
        waves = wave * A
        if clip:
            waves[waves > 10.] = 10.
            waves[waves < -10.] = -10.
        return(waves)



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
        Fs = self.PS.out_sampleFreq  # sample frequency
        stim = self.CPars['Stimulus']['Protocol']
#        level = None  # level is dbspl normally for models, but set to None for TDT (5V tone reference)
        seed = 32767
        #print('stim: ', stim)
        wave = None
        self.stim_vary = None
        self.total_trials = 1000
        if stim in ['Clicks']:
           wave = sound.ClickTrain(rate=Fs, duration=self.CPars['Stimulus']['Duration'], dbspl=level,
                            click_duration=self.CPars['Clicks']['Duration'], 
                            click_starts=1e-3*np.arange(self.CPars['Stimulus']['Delay']*1000.,
                                self.CPars['Clicks']['Interval']*self.CPars['Clicks']['Number']
                                    +self.CPars['Stimulus']['Delay'],
                                self.CPars['Clicks']['Interval']))

        elif stim in ['Tone RI', 'Tone Search']:
            if freq is None:
                freq = self.CPars['Stimulus']['Tone Frequency']*1000.
            if stim in ['Tone RI']:
                print (self.CPars['Stimulus'].keys())
                self.stim_vary = {'Intensity': Utility.seqparse(self.CPars['Stimulus']['Intensities'])[0][0]}
                self.total_trials = len(self.stim_vary['Intensity'])
            wave = sound.TonePip(rate=Fs, duration=self.CPars['Stimulus']['Duration']+self.CPars['Stimulus']['Delay'],
                            f0=freq, dbspl=level, 
                            pip_duration=self.CPars['Stimulus']['Duration'], pip_start=[self.CPars['Stimulus']['Delay']],
                            ramp_duration=self.CPars['Stimulus']['Rise-Fall']/1000)

        elif stim in ['One Tone']:
            self.total_trials = 1
            wave = sound.TonePip(rate=Fs, duration=self.CPars['Stimulus']['Duration']+self.CPars['Stimulus']['Delay'],
                            f0=self.tone_frequency*1000, dbspl=level, 
                            pip_duration=self.CPars['Stimulus']['Duration'], pip_start=[self.CPars['Stimulus']['Delay']],
                            ramp_duration=self.CPars['Stimulus']['Rise-Fall']/1000)

        elif stim in ['Tone SAM']:
            if freq is None:
                freq = self.CPars['Stimulus']['Tone Frequency']*1000.
            wave = sound.SAMTone(rate=Fs, duration=self.CPars['Stimulus']['Duration']+self.CPars['Stimulus']['Delay'],
                            f0=freq, dbspl=level, 
                            pip_duration=self.CPars['Stimulus']['Duration'], pip_start=[self.CPars['Stimulus']['Delay']],
                            ramp_duration=self.CPars['Stimulus']['Rise-Fall']/1000.,
                            fmod=self.CPars['Modulation/CMMR']['Frequency'],
                            dmod=self.CPars['Modulation/CMMR']['Depth'], seed=seed)

        elif stim in ['FM Sweep']:
            wave = sound.FMSweep(rate=Fs, duration=self.CPars['FMSweep']['Duration'], dbspl=level,
                                start=self.CPars['Stimulus']['Delay'], ramp=self.CPars['FMSweep']['Ramp Type'],
                                freqs=[self.CPars['FMSweep']['Freq Start']*1000.,
                                self.CPars['FMSweep']['Freq End']*1000.])

        elif stim in ['Noise RI', 'Noise Search']:
            if stim in ['Noise RI']:
                self.stim_vary = {'Intensity': Utility.seqparse(self.CPars['Stimulus']['Intensities'])[0][0]}
                self.total_trials = len(self.stim_vary['Intensity'])
            wave = sound.NoisePip(rate=Fs, duration=self.CPars['Stimulus']['Duration']+self.CPars['Stimulus']['Delay'],
                            f0=self.CPars['Stimulus']['Tone Frequency']*1000., dbspl=level, 
                            pip_duration=self.CPars['Stimulus']['Duration'], pip_start=[self.CPars['Stimulus']['Delay']],
                            ramp_duration=self.CPars['Stimulus']['Rise-Fall']/1000.,
                            fmod=self.CPars['Modulation/CMMR']['Frequency'], dmod=0., seed=seed)           
        elif stim in ['Noise SAM']:
            wave = sound.SAMNoise(rate=Fs, duration=self.CPars['Stimulus']['Duration']+self.CPars['Stimulus']['Delay'],
                            f0=self.CPars['Stimulus']['Tone Frequency']*1000., dbspl=level, 
                            pip_duration=self.CPars['Stimulus']['Duration'], pip_start=[self.CPars['Stimulus']['Delay']],
                            ramp_duration=self.CPars['Stimulus']['Rise-Fall']/1000.,
                            fmod=self.CPars['Modulation/CMMR']['Frequency'],
                            dmod=self.CPars['Modulation/CMMR']['Depth'], seed=seed)

        elif stim in ['FRA']: # frequency response area
            splseq = Utility.seqparse(self.CPars['Stimulus']['Intensities'])[0][0]
            freqseq = Utility.seqparse(self.CPars['Stimulus']['Frequencies'])[0][0]
            mat_spl, mat_freq = np.meshgrid(splseq, freqseq)
            self.stim_vary = {'Intensity': mat_spl.ravel(), 'Frequency': mat_freq.ravel()}
            self.total_trials = len(mat_spl.ravel())
            self.last_freq = None
            # note that for this one, the tone is computed at time of use

        elif stim in ['Noise Bands']:
            self.stim_vary = {'Intensity': Utility.seqparse(self.CPars['Stimulus']['Intensities'])[0][0]}
            self.total_trials = len(self.stim_vary['Intensity'])
            wave = sound.NoiseBandPip(rate=Fs,
                            duration=self.CPars['Stimulus']['Duration']+self.CPars['Stimulus']['Delay'],
                            f0=self.CPars['Stimulus']['Tone Frequency']*1000., dbspl=level, 
                            pip_duration=self.CPars['Stimulus']['Duration'], pip_start=[self.CPars['Stimulus']['Delay']],
                            ramp_duration=self.CPars['Stimulus']['Rise-Fall']/1000.,
                            seed=seed,
                            type=self.CPars['Noise Bands']['Type'],
                            noisebw=self.CPars['Noise Bands']['Noise BW']*1000.,
                            notchbw=self.CPars['Noise Bands']['Notch BW']*1000.,
                            centerfreq=self.CPars['Noise Bands']['CF']*1000.,
                            )
            
        elif stim in ['DMR']:
            wave = sound.DynamicRipple(rate=Fs, duration=5.0)
        
        elif stim in ['CMMR']:  # flanking type is "noise" (modulated), or "MultiTone", or "None".
                                # flankingPhase is comodulated or codeviant or random (if type is not None)
                                # spacing is band spacing in octaves (for flanking bands)
                                # 
            wave = sound.ComodulationMasking(rate=Fs, duration=self.CPars['Stimulus']['Duration']+
                                            self.CPars['Stimulus']['Delay'],
                        pip_duration=self.CPars['Stimulus']['Duration'],
                        pip_start=[self.CPars['Stimulus']['Delay']],
                        f0=self.CPars['Stimulus']['Tone Frequency']*1000.,
                        ramp_duration=self.CPars['Stimulus']['Rise-Fall']/1000.,
                        dbspl=level,
                        fmod=self.CPars['Modulation/CMMR']['Frequency'],
                        dmod=self.CPars['Modulation/CMMR']['Depth'],
                        flanking_type=self.CPars['Modulation/CMMR']['CMMR Flanking Type'],
                        flanking_spacing=self.CPars['Modulation/CMMR']['CMMR Flanking Spacing'],
                        flanking_phase=self.CPars['Modulation/CMMR']['CMMR Flanking Phase'],
                        flanking_bands=self.CPars['Modulation/CMMR']['CMMR Flanking Bands'],
                )

        elif stim in ['SSN']: # speech shaped noise
            # read the file:
            fname = "testsentence.wav"
            (rate, sig) = wav.read(fname) 
            duration = float(sig.shape[0] - 1)/rate
            wave = sound.SpeechShapedNoise(rate=Fs, duration=duration, waveform=sig[:,0], samplingrate=rate)

        elif stim in ['RSS']:
            wave = sound.RandomSpectrumShape(rate=Fs, duration=0.5, dbspl=level,
                        ramp='linear', ramp_duration=1e-2, f0=self.CPars['RSS Params']['CF']*1000,
                        pip_duration=self.CPars['Stimulus']['Duration'],
                        pip_start=[self.CPars['Stimulus']['Delay']],
                        amp_group_size=self.CPars['RSS Params']['Grouping'],
                        amp_sd=self.CPars['RSS Params']['Level SD'],
                        spacing=self.CPars['RSS Params']['Spacing'],
                        octaves=self.CPars['RSS Params']['Octaves'])  
        
        if stim in ['Noise Search', 'Tone Search']:
            self.searchmode = True  # search mode just runs "forever", until a stop is hit
        else:
            self.searchmode = False

        if wave is not None:
            self.wavesound = wave
            self.wave = self.map_voltage(stim, self.wavesound.sound, clip=True) # force computation and rescale and clip the waveform

    def show_wave(self):
        """
        Plot the waveform in the top graph
        """
        self.clearErrMsg()
        self.prepare_run()  # force computation/setup of stimulus
        self.plots['Wave'].clear()
        self.plots['Wave'].plot(self.wavesound.time, self.wave)
    
    def show_spectrogram(self):
        """
        Plot the spectrum in the middle graph
        
        If spectimage is checked in the main gui, also plot the spectrogram
        in a matplotlib window (must be closed to continue; is blocking)
        """
        self.clearErrMsg()
        self.prepare_run()
        Fs = self.PS.out_sampleFreq
        # show the long term spectrum.
        f, Pxx_spec = scipy.signal.periodogram(self.wave, Fs) #, window='flattop', nperseg=8192,
                       # noverlap=512, scaling='spectrum')
        self.plots['LongTermSpec'].clear()
        self.plots['LongTermSpec'].plot(f[1:], np.sqrt(Pxx_spec)[1:], pen=pg.mkPen('y'))
        #self.plots['LongTermSpec'].setLogMode(x=True, y=False)

        print (self.maingui.spectimage)
        if self.maingui.spectimage:  # enable spectrogram plot
            import matplotlib.pyplot as mpl
            ax1 = mpl.subplot(211)
            mpl.plot(self.wavesound.time, self.wave)
            mpl.subplot(212, sharex=ax1)
            Pxx, freqs, bins, im = mpl.specgram(self.wave, NFFT=128, Fs=Fs, noverlap=64, pad_to=256)
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
            #self.img.setImage(Sxx.T)
    
    def clearErrMsg(self):
        """
        Reset the error notificatoin to the standard ready indicator
        """
        self.maingui.permStatusMessage.setText('<b><font color="#00FF00">Ready</b>')


# Build GUI and window

class BuildGui():
    def __init__(self):
        self.app = pg.mkQApp()
        self.mainwin = QtGui.QMainWindow()
        self.win = QtGui.QWidget()
        self.layout = QtGui.QGridLayout()
        self.win.setLayout(self.layout)
        self.mainwin.setCentralWidget(self.win)
        self.mainwin.show()
        self.mainwin.setWindowTitle('Stim Controller')
        self.mainwin.setGeometry( 100 , 100 , 1024 , 800)
        self.spectimage = False
        self.TDTTankDirectory = ''
        self.TT = TDT.TDTTankInterface()
        print('self.TT.available: ', self.TT.available)
        self.statusBar = QtGui.QStatusBar()
        self.mainwin.setStatusBar(self.statusBar)
        self.statusMessage = QtGui.QLabel('')
        self.statusBar.addWidget(self.statusMessage)
        self.permStatusMessage = QtGui.QLabel('<b><font color="#00FF00">Ready</b>')
        self.statusBar.addPermanentWidget(self.permStatusMessage)
        

        # retrieve recent path
        self.configfilename = 'config.ini'
        if not os.path.isfile(self.configfilename):
            # create a configuration file
            parser = ConfigParser.SafeConfigParser()
            # initialize parser
            parser.add_section('TDTTanks')
            parser.set('TDTTanks', 'dir', '')
            fh = open(self.configfilename, 'w')
            parser.write(fh)
            fh.close()
        else:
            parser = ConfigParser.SafeConfigParser()
            parser.read('config.ini')
            self.TDTTankDirectory = parser.get('TDTTanks', 'dir')
            print('tankd dir: ', self.TDTTankDirectory)
            self.TT.tank_directory = self.TDTTankDirectory
            if len(self.TDTTankDirectory) > 0:
                self.TT.open_tank()
                lastblock = self.TT.find_last_block()
                self.TT.close_tank()
                self.TT.show_tank_path()

        # Define parameters that control aquisition and buttons...
        params = [
            {'name': 'Stimulus', 'type': 'group', 'children': [
                {'name': 'Protocol', 'type': 'list', 'values': ['Noise Search', 'Tone Search',
                        'Tone RI', 'Noise RI', 'FRA', 'Clicks', 
                        'CMMR', 'RSS', 'DMR', 'SSN',
                        'Tone SAM', 'Noise SAM', 'FM Sweep',
                        'Noise Bands'],
                        'value': 'Noise Search'},
                {'name': 'Tone Frequency', 'type': 'float', 'value': 4.0, 'step': 1.0, 'limits': [0.5, 99.0],
                    'suffix': 'kHz', 'default': 4.0},
                {'name': 'Attenuator', 'type': 'float', 'value': 50, 'step': 5.0, 'limits': [0., 120.0],
                    'suffix': 'dB', 'default': 50.0},
                {'name': 'Rise-Fall', 'type': 'float', 'value': 2.5, 'step': 0.5, 'limits': [0.5, 20.],
                    'suffix': 'ms', 'default': 2.5},
                {'name': 'Intensities', 'type': 'str', 'value': '90;20/-10',
                    'suffix': 'dBAttn', 'default': '90;20/-10'},
                {'name': 'Frequencies', 'type': 'str', 'value': '4;48/8l',
                    'suffix': 'kHz', 'default': '4;48/8l'},
                                            
                {'name': 'Repetitions', 'type': 'int', 'value': 1, 'limits': [1, 10000], 'default': 1,
                     'tip': 'Number of Stimuli per sweep'},
                {'name': 'Intertrial Interval', 'type': 'float', 'value': 1., 'limits': [0.5, 300.], 
                    'suffix': 's', 'default': 1.0, 'tip': 'Time between sweeps (trials) in FRA and RI protocols'},
                {'name': 'Interstimulus Interval', 'type': 'float', 'value': 0.8, 'limits': [0.01, 300.], 
                    'suffix': 's', 'default': 0.8, 'tip': 'Time between stimuli in a sweep'},
                {'name': 'Randomize', 'type': 'bool', 'value': False, 'default': False,
                    'tip': 'Randomize presentation order in all dimensions'},
                {'name': 'Duration', 'type': 'float', 'value': 0.2, 'step': 0.05, 'limits': [0.001, 10],
                    'suffix': 's', 'default': 0.2, 'tip': 'Sound duration, in seconds'},
                {'name': 'Delay', 'type': 'float', 'value': 0.01, 'step': 0.05, 'limits': [0.001, 10],
                    'suffix': 's', 'default': 0.2, 'tip': 'Sound delay, in seconds'},
             ]},
             {'name': 'Clicks', 'type': 'group', 'expanded': False, 'children': [
                  {'name': 'Interval', 'type': 'float', 'value': 50., 'step': 5.0, 'limits': [1., 1000.0],
                    'suffix': 'ms', 'default': 50.0},
                  {'name': 'Number', 'type': 'int', 'value': 4, 'step': 1, 'limits': [1, 200.0],
                    'default': 4},
                  {'name': 'Duration', 'type': 'float', 'value': 1e-4, 'step': 10e-6, 
                      'limits': [10e-6, 1e-3], 'suffix': 's', 'default': 1e-4},
             ]},

             {'name': 'FMSweep', 'type': 'group', 'expanded': False, 'children': [
                  {'name': 'Duration', 'type': 'float', 'value': 0.5, 'step': 0.05, 
                      'limits': [5e-3, 10], 'suffix': 's', 'default': 0.5},
                  {'name': 'Ramp Type', 'type': 'list', 'values': ['linear', 'logarithmic'], 'value': 'linear'},
                  {'name': 'Freq Start', 'type': 'float', 'value': 4, 'step': 1, 'limits': [1., 100.0],
                    'default': 4},
                  {'name': 'Freq End', 'type': 'float', 'value': 48, 'step': 1, 'limits': [1., 100.0],
                    'default': 48},
             ]},

             {'name': 'Modulation/CMMR', 'type': 'group', 'expanded': False, 'children': [
                  {'name': 'Frequency', 'type': 'float', 'value': 40., 'step': 5.0, 'limits': [1., 1000.0],
                    'suffix': 'Hz', 'default': 40.0},
                  {'name': 'Depth', 'type': 'float', 'value': 50.0, 'step': 5.0, 'limits': [0.0, 200.0],
                    'suffix': '%', 'default': 50.},
                  {'name': 'CMMR Flanking Type', 'type': 'list', 'values': ['None', 'MultiTone', 'NBnoise'], 'value': 'MultiTone'},
                  {'name': 'CMMR Flanking Phase', 'type': 'list', 'values': ['Comodulated', 'Codeviant', 'Random'],
                   'value': 'Comodulated'},
                  {'name': 'CMMR Flanking Bands', 'type': 'int', 'value': 2, 'step': 1, 'limits': [0, 10],
                   'default': 2},
                  {'name': 'CMMR Flanking Spacing', 'type': 'float', 'value': 1, 'step': 1/8., 
                      'limits': [1/16., 2.], 'suffix': 'octaves', 'default': 1},
             ]},
             {'name': 'RSS Params', 'type': 'group', 'expanded': False,  'children': [
                 {'name': 'CF', 'type': 'float', 'value': 16.0, 'step': 2.0, 'limits': [1.0, 64.],
                 'suffix': 'kHz', 'default': 16.0},
                 {'name': 'Grouping', 'type': 'int', 'value': 8, 'step': 1, 'limits': [1, 64],
                   'default': 8},
                 {'name': 'Level SD', 'type': 'float', 'value': 12.0, 'step': 2.0, 'limits': [0.0, 40.],
                 'suffix': 'dB', 'default': 12.0},
                 {'name': 'Spacing', 'type': 'int', 'value': 64, 'step': 2, 'limits': [1, 128],
                 'suffix': '/octave', 'default': 64},
                 {'name': 'Octaves', 'type': 'float', 'value': 3.0, 'step': 0.5, 'limits': [0.5, 16.],
                   'default': 3.0},
             ]},
             {'name': 'Noise Bands', 'type': 'group', 'expanded': False, 'children': [
                 {'name': 'Type', 'type': 'list', 'values': ['Bandpass', 'BP+Notch'],
                  'value': 'Bandpass'},
                 {'name': 'Notch BW', 'type': 'float', 'value': 1.0, 'step': 1.0, 'limits': [0.05, 10.],
                 'suffix': 'kHz', 'default': 1.0},
                 {'name': 'Noise BW', 'type': 'float', 'value': 10.0, 'step': 1.0, 'limits': [1.0, 64],
                 'suffix': 'kHz', 'default': 10.0},
                 {'name': 'CF', 'type': 'float', 'value': 5.0, 'limits': [0.1, 40],
                 'suffix': 'kHz', 'default': 5.0},
             ]},
                    

            {'name': 'File From Disk', 'type': 'str', 'value': 'test.wav', 'default': 'test.wav'},

        ]
        
        self.ptree = ParameterTree()
        self.ptreedata = Parameter.create(name='params', type='group', children=params)
        self.ptree.setParameters(self.ptreedata)
        ptreewidth = 120

        #print (dir(self.ptreedata))
        
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
        self.btn_waveform = QtGui.QPushButton("Wave")
        self.btn_spectrum = QtGui.QPushButton("Spectrum")
        self.btn_run = QtGui.QPushButton("Run")
        self.btn_pause = QtGui.QPushButton("Pause")
        self.btn_continue = QtGui.QPushButton("Continue")
        self.btn_stop = QtGui.QPushButton("Stop")
        self.btn_quit = QtGui.QPushButton("Quit")
        self.btn_tdt = QtGui.QPushButton("TDT Tank")
        self.label_status = QtGui.QLabel('Stopped')
        self.label_trialctr = QtGui.QLabel('Trial: 0')
        self.label_status.sizeHint = QtCore.QSize(100, 20)
        self.label_trialctr.setAutoFillBackground(True)
        self.label_trialctr.sizeHint = QtCore.QSize(100, 20)
        self.spect_check = QtGui.QCheckBox('Spectrogram')
        self.spect_check.setChecked(False)  # just be sure.
        hbox = QtGui.QGridLayout()
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
        self.layout.addWidget(self.ptree, 1, 0, 4, 2) # Parameter Tree on left
        self.layout.setColumnMinimumWidth(0, ptreewidth)
 
        # add space for the graphs
        view = pg.GraphicsView()
        glayout = pg.GraphicsLayout(border=(50,50,50))
        view.setCentralItem(glayout)
        self.layout.addWidget(view, 0, 2, 5, 3)  # data plots on right
        self.plots = {}
        self.plots['Wave'] = glayout.addPlot()
        self.plots['Wave'].getAxis('left').setLabel('V', color="#ff0000")
        self.plots['Wave'].setTitle('Waveform', color="#ff0000")
        self.plots['Wave'].getAxis('bottom').setLabel('t (sec)', color="#ff0000")
        self.plots['Wave'].setYRange(-1, 1)

        glayout.nextRow()
        self.plots['LongTermSpec'] = glayout.addPlot()
        self.plots['LongTermSpec'].getAxis('left').setLabel('V', color="#ff0000")
        self.plots['LongTermSpec'].setTitle('LongTerm Spectrum', color="#ff0000")
        self.plots['LongTermSpec'].getAxis('bottom').setLabel('F (Hz)', color="#ff0000")


    #     if self.spectimage:
    #         self.img = pg.ImageView() # view=self.plots['Spec'])
    #         arr = np.random.random((100, 32))
    #         self.img.setImage(arr)
    #         self.img.ui.roiBtn.hide()
    # #        self.img.ui.menuBtn.hide()
    #         self.img.show()
    #     else:
        self.img = None

        glayout.nextRow()
        l2 = glayout.addLayout(colspan=3, border=(50,0,0))  # embed a new layout
        l2.setContentsMargins(10,10,10,10)
        self.plots['Plot1'] = l2.addPlot(Title="Plot1")
#        self.l2.addWidget(self.plots['Plot1'])
        self.plots['Plot1'].getAxis('bottom').setLabel('F (kHz)')
        self.plots['Plot1'].getAxis('left').setLabel('dB ATTN')
        self.plots['Plot1'].setTitle('FRA')
        self.plots['Plot1'].setXRange(0, 50, padding=0)
        #self.plots['Plot1'].setLogMode(x=True)
        self.plots['Plot1'].setYRange(125, -5, padding=0)
        xd = np.arange(2, 48, 1)
       # xd = np.logspace(np.log2(2), np.log2(64), 50, base=2)
       # print ('xd: ', xd)
        yd = np.arange(120, 5, -5)
        spots = []
        self.lastPoint = None
        for i in range(xd.shape[0]):
            for j in range(yd.shape[0]):
                spots.append({'pos': (xd[i], yd[j]), 'size': 7, 'pen': {'color': 'k', 'width': 0.5, 'alpha': 0.5},
                    'brush': pg.mkBrush('b')})
        self.spi = pg.ScatterPlotItem(size=7, pen=pg.mkPen('k'), brush=pg.mkBrush('b'), symbol='s')
        self.spi.addPoints(spots)
        self.plots['Plot1'].addItem(self.spi)
        self.spi.getViewBox().invertY(True)
        self.spi.sigClicked.connect(self.getClickedLocation)
        #cross hair
        # vLine = pg.InfiniteLine(angle=90, movable=True)
        # hLine = pg.InfiniteLine(angle=0, movable=True)
        # self.plots['Plot1'].addItem(vLine, ignoreBounds=False)
        # self.plots['Plot1'].addItem(hLine, ignoreBounds=False)
        # vb = self.plots['Plot1'].vb

        # def mouseMoved(evt):
        #     pos = evt[0]  ## using signal proxy turns original arguments into a tuple
        #     if self.plots['Plot1'].sceneBoundingRect().contains(pos):
        #         mousePoint = vb.mapSceneToView(pos)
        #         index = int(mousePoint.x())
        #         if index > 0 and index < len(data1):
        #             label.setText("<span style='font-size: 12pt'>x=%0.1f,   <span style='color: red'>y1=%0.1f</span>,   <span style='color: green'>y2=%0.1f</span>" % (mousePoint.x(), data1[index], data2[index]))
        #         vLine.setPos(mousePoint.x())
        #         hLine.setPos(mousePoint.y())
        # proxy = pg.SignalProxy(self.plots['Plot1'].scene().sigMouseMoved, rateLimit=60, slot=mouseMoved)   

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
        self.controller = Controller(self.ptreedata, self.plots, self.img, self)  # we pass the gui also

        self.controller.setAllParameters(params)  # synchronize parameters with the tree
        self.ptreedata.sigTreeStateChanged.connect(self.controller.change)  # connect parameters to their updates
        
        # now connect the buttons
        self.recentpath = ''
        self.btn_waveform.clicked.connect(self.controller.show_wave)
        self.btn_spectrum.clicked.connect(self.controller.show_spectrogram)
        self.btn_tdt.clicked.connect(self.getTDTTank) # self.TT.set_tank_path)
     #    self.ButtonEvents = QtCore.QTimer() # get a Q timer
     #    self.btn_stop.clicked.connect(timer, SIGNAL(timeout()), this, SLOT(processOneThing()));
     # timer->start();
     
        self.btn_run.clicked.connect(self.controller.start_run)


        self.btn_pause.clicked.connect(self.controller.pause_run)
        self.btn_continue.clicked.connect(self.controller.next_stimulus)
        self.btn_stop.clicked.connect(self.controller.stop_run)
        self.btn_quit.clicked.connect(self.controller.quit)
        self.spect_check.clicked.connect(self.speccheck)
        self.updateStatusMessage()

    def speccheck(self):
        self.spectimage = self.spect_check.isChecked()

    def getTDTTank(self, dirname=None):
        filedialog = QtGui.QFileDialog()
        filedialog.setFileMode(QtGui.QFileDialog.Directory)
        self.TT.tank_directory = str(filedialog.getExistingDirectory(None, "Select Tank Directory", self.recentpath,
                                    QtGui.QFileDialog.ShowDirsOnly))
        self.recentpath = self.TT.tank_directory
#        print('Tank dir selected: ', self.TT.tank_directory)
        self.setTankIni(self.TT.tank_directory)
        self.TT.open_tank()
        lastblock = self.TT.find_last_block()
        self.TT.close_tank()
        self.TT.show_tank_path()
        self.updateStatusMessage()

    def setTankIni(self, newtankname):
        parser = ConfigParser.SafeConfigParser()
        parser.read('config.ini')
        parser.set('TDTTanks', 'dir', newtankname)
        fh = open(self.configfilename, 'w')
        parser.write(fh)
        fh.close()

    def updateStatusMessage(self):
        if self.TT.available is False:
            message = ('No TDT Tank')
        else:
            message = ('Tank: {0:s}  CurrentBlock: {1:d}'.format(self.TT.tank_directory, self.TT.lastblock))
        self.statusMessage.setText(message)
        
    def getClickedLocation(self, points):
        # print (dir(points))
        # print (points.event())
        # print('mouse click: ', points.mouseClickEvent(points.event()))
        # print('mouse doubleclick: ', points.mouseDoubleClickEvent(points.event()))
        self.mousePoint = points.ptsClicked[0].viewPos()
        print('mousepoint: ', self.mousePoint.x(), self.mousePoint.y()) 
        points.ptsClicked[0].setBrush(pg.mkBrush('r'))
        points.ptsClicked[0].setSize(7)
        if self.lastPoint is None:
            self.lastPoint = points.ptsClicked[0]
        else:
            self.lastPoint.setBrush(pg.mkBrush('b'))
            self.lastPoint.setSize(7)
            self.lastPoint = points.ptsClicked[0]
        
        stimpars = self.ptreedata.param('Stimulus').items.keys()[0]  # force to One Tone mode
        stimpars.param.names['Protocol'].setValue('One Tone')
#        stimpars.param.emitStateChanged()  # trigger

        self.controller.protocol = stimpars.param.names['Protocol'].value()
        self.controller.tone_frequency = self.mousePoint.x()
        self.controller.attn = self.mousePoint.y()
       # self.controller.prepare_run()
        self.controller.start_run()
    
if __name__ == '__main__':
    gui = BuildGui()
    
    ## Start Qt event loop unless running in interactive mode.
    ## Event loop will wait for the GUI to activate the updater and start sampling.
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
    gui.controller.quit()