from __future__ import print_function

import sys
import datetime
import numpy as np
import scipy.signal
import scipy.io.wavfile as wav
from collections import OrderedDict
import pyqtgraph as pg
from PyQt4 import QtGui, QtCore
from pyqtgraph.parametertree import Parameter, ParameterTree
import pystim
import Utility
import sound
import pprint

pp = pprint.PrettyPrinter(indent=4)


class Controller(object):
    def __init__(self, ptreedata, plots, img):
        self.Vscale = 5.0  # CRITICAL: voltate for staqndard tone level
        self.PS = pystim.PyStim(hdw=['PA5', 'NIDAQ', 'RZ5D'])
        self.ptreedata = ptreedata
        self.plots = plots  # access to plotting area
        self.img = img
         
        
        self.CPars = OrderedDict()
        
#        self.setAllParameters(ptreedata)
        
        self.CPars['tone_frequency'] = 4.0  # tone pip frequency, Hz
        self.CPars['duration'] = 0.2  # stimulus duration, ms
        self.CPars['delay'] = 0.01  # delay to start of stimulus, s
        self.CPars['RI'] = '20;100/10'  # seqparse for Rate-Intensity intensity series
        self.CPars['frequencies'] = '4;48/8l'  # seqparse for tone pip FRA
        self.CPars['dMod'] = 50.  # modulation depth
        self.CPars['fMod'] = 40. # modulation frequency
        self.CPars['RF'] = 2.5  # rise-fall time, ms
        self.CPars['nreps'] = 5  # repetitions in a single sweep
        self.CPars['interstimulus_interval'] = 1.0  # time between stimuli, s
        self.CPars['protocol'] = 'Tone RI'
        self.CPars['randomize'] = False  # randomize order of presentation (or not...)
        
        self.CPars['click_interval'] = 50.
        self.CPars['click_N'] = 4
        self.CPars['click_duration'] = 1e-4
        
        self.CPars['RSS_cf'] = 16.0
        self.CPars['RSS_grouping'] = 8
        self.CPars['RSS_octaves'] = 3
        self.CPars['RSS_SD'] = 12.
        self.CPars['RSS_spacing'] = 64
        
        self.CPars['CMMR_flanking_bands'] = 3  # flanking bands
        self.CPars['CMMR_flanking_phase'] = 'comodulated'  # flanking bands comodulated 
        self.CPars['CMMR_flanking_spacing'] = 0.5  # octaves
        
                
    # def setAllParameters(self, params):
    #     """
    #     Set all of the local parameters from the parameter tree
    #
    #     Parameters
    #     ----------
    #     ptree : ParameterTree object
    #
    #     Returns
    #     -------
    #     Nothing
    #     """
    #     pass
    

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
#        print('changes: ', changes)
        for param, change, data in changes:
            path = self.ptreedata.childPath(param)
#            print( 'path: ', path)
            # if path is not None:
            #     childName = '.'.join(path)
            # else:
            #     childName = param.name()
            #
            # Parameters and user-supplied information
            #
            #print('Path[1]: ', path[1])

            if path[0] == 'Stimulus':
                if path[1] == 'Protocol':
                    self.CPars['protocol'] = data
                if path[1] == 'Repetitions':
                    self.CPars['nreps'] = data
                if path[1] == 'Interstimulus Interval':
                    self.CPars['interstimulus_interval'] = data
                if path[1] == 'Duration':
                    self.CPars['duration'] = data
                if path[1] == 'Delay':
                    self.CPars['delay'] = data
                if path[1] == 'Tone Frequency':
                    self.CPars['tone_frequency'] = data
                if path[1] == 'Rise-Fall':
                    self.CPars['RF'] = data
                if path[1] == 'Intensities':
                    self.CPars['RI'] = data
                if path[1] == 'Frequencies':
                    self.CPars['frequencies'] = data

            if path[0] == 'Clicks':
                if path[1] == 'Interval':
                    self.CPars['click_interval'] = data
                if path[1] == 'Number':
                    self.CPars['click_N'] = data
                if path[1] == 'Duration':
                    self.CPars['click_duration'] = data
            
            if path[0] == 'Modulation':
                if path[1] == 'Modulation Depth':
                    self.CPars['dMod'] = data
                if path[1] == 'Modulation Frequency':
                    self.CPars['fMod'] = data
                if path[1] == 'CMMR Flanking Type':
                    self.CMMR_flanking_type = data
                if path[1] == 'CMMR Flanking Bands':
                    self.CPars['CMMR_flanking_bands'] = data
                if path[1] == 'CMMR Flanking Phase':
                    self.CPars['CMMR_flanking_phase'] = data
                if path[1] == 'CMMR Flanking Spacing':
                    self.CPars['CMMR_flanking_spacing'] = data
                    
            if path[0] == 'RSS Params':
                if path[1] == 'CF':
                    self.CPars['RSS_cf'] = data
                if path[1] == 'Grouping':
                    self.CPars['RSS_grouping'] = data
                if path[1] == 'Octaves':
                    self.CPars['RSS_octaves'] = data
                if path[1] == 'Level SD':
                    self.CPars['RSS_SD'] = data
                if path[1] == 'Spacing':
                    self.CPars['RSS_spacing'] = data
    
    def showParameters(self):
        pp.pprint(self.CPars)
    
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
        for p in params[0]['children']:
            if p['type'] == 'action':
                continue
            if p == 'Repetitions':
                self.CPars['nreps'] = p['value']
            if p == 'Interstimulus Interval':
                self.CPars['interstimulus_interval'] = p['value']
            if p == 'Duration':
                self.CPars['duration'] = p['value']
            if p == 'Delay':
                self.CPars['delay'] = p['value']
            if p == 'Tone Frequency':
                self.CPars['tone_frequency'] = p['value']
            if p == 'Rise-Fall':
                self.CPars['RF'] = p['value']
            if p == 'Intensities':
                self.CPars['RI'] = p['value']
            if p == 'Frequencies':
                self.CPars['frequencies'] = p['value']
            if p == 'Modulation Depth':
                self.CPars['dMod'] = p['value']
            if p == 'Modulation Frequency':
                self.CPars['fMod'] = p['value']

    def run(self):
        """
        Initialize variables for the start of a run
        then use continueRun to get timers
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Nothing
        """

        self.showParameters()
        
        self.runtime = 0
        self.NSamples = 0
        self.running = True
        self.stop_hit = False
        self.startTime = datetime.datetime.now()
        self.prepare_run()  # reset the data arrays
        self.continue_run()

    def pause_run(self):
        self.pause_hit = True
        
    def continue_run(self):
        """
        Start the timing and data collection, without reinitializing
        any variables
        Parameters
        ----------
        None
        
        Returns
        -------
        Nothing
        
        """
        self.timer = QtCore.QTimer()
        # self.timer.timeout.connect(updater.update)
        self.timedWrite = QtCore.QTimer()
        # self.timedWrite.timeout.connect(updater.storeData)
        # self.update() # do the first update, then start time
        # self.timer.start(self.readInterval * 1000)
        # self.timedWrite.start(30 * 1000.)  # update file in 1 minute increments
        # print ('Playing %s' % self.CPars['protocol'])
        # do diferently according to protocol:
        
        protocol = self.CPars['protocol']
        if protocol in ['Noise Search', 'Tone Search']:
            spl=80.
            self.prepare_run(level=spl)  # same stimulus all the time, just compute
            while not self.stop_hit:
                self.PS.play_sound(self.wave.sound*self.Vscale, self.wave.sound*self.Vscale,
                    samplefreq=self.PS.out_sampleFreq,
                    isi=self.CPars['interstimulus_interval'], reps=self.CPars['nreps'], attns=self.convert_spl_attn(spl))
        elif protocol in ['Tone RI', 'Noise RI']:
            seq = Utility.seqparse(self.CPars['RI'])
            splseq = seq[0][0]
            for spl in splseq:
                self.prepare_run(level = spl)  # same stimulus all the time, just amplitude changes
                print('spl:', spl)
                print('Protocol {0:s}  attn: {1:3.1f}'.format(protocol, spl))
                self.PS.play_sound(self.wave.sound*self.Vscale, self.wave.sound*self.Vscale,
                    samplefreq=self.PS.out_sampleFreq,
                    isi=self.CPars['interstimulus_interval'], reps=self.CPars['nreps'], attns=self.convert_spl_attn(spl))
                if self.stop_hit:
                    break

        elif protocol in ['FRA']:
            splseq = Utility.seqparse(self.CPars['RI'])[0][0]
            freqseq = Utility.seqparse(self.CPars['frequencies'])[0][0]
            print('spls: ', self.CPars['RI'], splseq)
            print('freq: ', self.CPars['frequencies'], freqseq)
            for spl in splseq:
                for freq in freqseq:
                    self.prepare_run(freq=freq*1000., level=spl)  # same stimulus all the time, just amplitude changes
                    print('Protocol {0:s}  freq: {1:6.3f}  spl: {2:3.1f}'.format(protocol, freq, spl))
                    self.PS.play_sound(self.wave.sound*self.Vscale, self.wave.sound*self.Vscale,
                        samplefreq=self.PS.out_sampleFreq,
                        isi=self.CPars['interstimulus_interval'], reps=self.CPars['nreps'],
                        attns=self.convert_spl_attn(spl))
                    if self.stop_hit:
                        break

        else:
            spl = 80
            self.prepare_run(level = spl)  # same stimulus all the time, just amplitude changes
            self.PS.play_sound(self.wave.sound*self.Vscale, self.wave.sound*self.Vscale,
                samplefreq=self.PS.out_sampleFreq,
                isi=self.CPars['interstimulus_interval'], reps=self.CPars['nreps'], attns=self.convert_spl_attn(spl))
            
        
        
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
        self.timer.stop()
        self.stop_hit = True
        # self.timedWrite.stop()
#        self.storeData()
    
    
    def quit(self):
        self.PS.HwOff()
        exit(0)
    
    def convert_spl_attn(self, spl):
        return [100.-spl, 100.-spl]  # just crude, need to clean up

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
        stim = self.CPars['protocol']
#        level = None  # level is dbspl normally for models, but set to None for TDT (5V tone reference)
        seed = 32767
        #print('stim: ', stim)
        wave = None
        if stim in ['Clicks']:
           wave = sound.ClickTrain(rate=Fs, duration=self.CPars['duration'], dbspl=level,
                            click_duration=self.CPars['click_duration'], 
                            click_starts=1e-3*np.linspace(self.CPars['delay'], 
                            self.CPars['click_interval']*self.CPars['click_N'], self.CPars['click_interval']))

        elif stim in ['Tone RI', 'Tone Search', 'FRA']:
            if freq is None:
                freq = self.CPars['tone_frequency']*1000.
            wave = sound.TonePip(rate=Fs, duration=self.CPars['duration']+self.CPars['delay'],
                            f0=freq, dbspl=level, 
                            pip_duration=self.CPars['duration'], pip_start=[self.CPars['delay']],
                            ramp_duration=self.CPars['RF']/1000)
        elif stim in ['Tone SAM']:
            if freq is None:
                freq = self.CPars['tone_frequency']*1000.
            wave = sound.SAMTone(rate=Fs, duration=self.CPars['duration']+self.CPars['delay'],
                            f0=freq, dbspl=level, 
                            pip_duration=self.CPars['duration'], pip_start=[self.CPars['delay']],
                            ramp_duration=self.CPars['RF']/1000.,
                            fmod=self.CPars['fMod'], dmod=self.CPars['dMod'], seed=seed)

        elif stim in ['FM Sweep']:
            wave = sound.FMSweep(rate=Fs, duration=0.5, dbspl=level,
                               start=0., ramp='linear', freqs=[16000, 200])

        elif stim in ['Noise RI', 'Noise Search']:
            wave = sound.NoisePip(rate=Fs, duration=self.CPars['duration']+self.CPars['delay'],
                            f0=self.CPars['tone_frequency']*1000., dbspl=level, 
                            pip_duration=self.CPars['duration'], pip_start=[self.CPars['delay']],
                            ramp_duration=self.CPars['RF']/1000.,
                            fmod=self.CPars['fMod'], dmod=0., seed=seed)           
        elif stim in ['Noise SAM']:
            wave = sound.SAMNoise(rate=Fs, duration=self.CPars['duration']+self.CPars['delay'],
                            f0=self.CPars['tone_frequency']*1000., dbspl=level, 
                            pip_duration=self.CPars['duration'], pip_start=[self.CPars['delay']],
                            ramp_duration=self.CPars['RF']/1000.,
                            fmod=self.CPars['fMod'], dmod=self.CPars['dMod'], seed=seed)

        elif stim in ['DMR']:
            wave = sound.DynamicRipple(rate=Fs, duration=5.0)
        
        elif stim in ['CMMR']:  # flanking type is "noise" (modulated), or "3Tone", or "None".
                                # flankingPhase is comodulated or codeviant or random (if type is not None)
                                # spacing is band spacing in octaves (for flanking bands)
                                # 
            wave = sound.ComodulationMasking(rate=Fs, duration=self.CPars['duration'],
                f0=self.CPars['tone_frequency']*1000., 
                dbspl=level, fmod=self.CPars['fMod'], dmod=self.CPars['dMod'],
                flanking_type=self.CMMR_flanking_type, flanking_spacing=self.CPars['CMMR_flanking_spacing'],
                flanking_phase=self.CPars['CMMR_flanking_phase'], flanking_bands=self.CPars['CMMR_flanking_bands'],
                )

        elif stim in ['SSN']: # speech shaped noise
            # read the file:
            fname = "testsentence.wav"
            (rate, sig) = wav.read(fname) 
            duration = float(sig.shape[0] - 1)/rate
            wave = sound.SpeechShapedNoise(rate=Fs, duration=duration, waveform=sig[:,0], samplingrate=rate)

        elif stim in ['RSS']:
            wave = sound.RandomSpectrumShape(rate=Fs, duration=0.5, dbspl=level,
                ramp='linear', ramp_duration=1e-2, f0=self.CPars['RSS_cf']*1000, pip_duration=0.4,
                pip_start=[self.CPars['delay']], amp_group_size=self.CPars['RSS_grouping'], amp_sd=self.CPars['RSS_SD'],
                spacing=self.CPars['RSS_spacing'], octaves=self.CPars['RSS_octaves'])  

        self.wave = wave # rescale the waveform

    def show_wave(self):
        self.prepare_run()  # force computation/setup of stimulus
        self.plots['Wave'].clear()
        self.plots['Wave'].plot(self.wave.time, self.wave.sound*self.Vscale)
    
    # redundant: pyqtgraph provides spectrum plot
    # def show_spectrum(self):
    #     self.prepare_run()
    #     Fs = self.PS.out_sampleFreq
    #     f, Pxx_spec = scipy.signal.periodogram(self.wave.sound, Fs) #, window='flattop', nperseg=8192,
    #                        # noverlap=512, scaling='spectrum')
    #     self.plots['Wave'].plot(f[1:], np.sqrt(Pxx_spec)[1:])
    #     self.plots['Wave'].setLogMode(x=True, y=True)

    def show_spectrogram(self):
        self.prepare_run()
        Fs = self.PS.out_sampleFreq
        specfreqs, spectime, Sxx = scipy.signal.spectrogram(self.wave.sound*self.Vscale, nperseg=int(0.01*Fs), fs=Fs)
        thr = 0. # 1e-8
        Sxx[Sxx <= thr] = thr
        # pos = np.array([0., 1., 0.5, 0.25, 0.75])
        # color = np.array([[0,255,255,255], [255,255,0,255], [0,0,0,255], (0, 0, 255, 255), (255, 0, 0, 255)], dtype=np.ubyte)
        # cmap = pg.ColorMap(pos, color)
        # lut = cmap.getLookupTable(0.0, 1.0, 256)
        # # set colormap
        # # print (dir(self.img))
        # # print (dir(self.img.imageItem))
        # self.img.imageItem.setLookupTable(lut)
#        self.img.setLevels([-40,50]) 
        self.img.setImage(Sxx.T)
        
        # also show the long term spectrum.
        f, Pxx_spec = scipy.signal.periodogram(self.wave.sound*self.Vscale, Fs) #, window='flattop', nperseg=8192,
                       # noverlap=512, scaling='spectrum')
        self.plots['LongTermSpec'].clear()
        self.plots['LongTermSpec'].plot(f[1:], np.sqrt(Pxx_spec)[1:], pen=pg.mkPen('y'))
        #self.plots['LongTermSpec'].setLogMode(x=True, y=False)
        

        
# Build GUI and window

class BuildGui(object):
    def __init__(self):
        self.app = pg.mkQApp()
        self.win = QtGui.QWidget()
        self.layout = QtGui.QGridLayout()
        self.win.setLayout(self.layout)
        self.win.show()
        self.win.setWindowTitle('Stim Controller')
        self.win.setGeometry( 100 , 100 , 1024 , 800)

        # Define parameters that control aquisition and buttons...
        params = [
            {'name': 'Stimulus', 'type': 'group', 'children': [
                {'name': 'Protocol', 'type': 'list', 'values': ['Noise Search', 'Tone Search',
                        'Tone RI', 'Noise RI', 'FRA',
                        'RSS', 'DMR', 'SSN', 'Tone SAM', 'Noise SAM', 'Clicks', 'FM Sweep',
                        'NotchNoise', 'BandPass Noise'], 'value': 'Tone RI'},
                {'name': 'Tone Frequency', 'type': 'float', 'value': 4.0, 'step': 1.0, 'limits': [0.5, 99.0],
                    'suffix': 'kHz', 'default': 4.0},
                {'name': 'Rise-Fall', 'type': 'float', 'value': 2.5, 'step': 0.5, 'limits': [0.5, 20.],
                    'suffix': 'ms', 'default': 2.5},
                {'name': 'Intensities', 'type': 'str', 'value': '20;90/10',
                    'suffix': 'dB', 'default': '20;90/10'},
                {'name': 'Frequencies', 'type': 'str', 'value': '4;48/8l',
                    'suffix': 'kHz', 'default': '4;48/8l'},
                                            
                {'name': 'Repetitions', 'type': 'int', 'value': 5, 'limits': [1, 10000], 'default': 5, 'tip': 'Stimuli per sweep'},
                {'name': 'Interstimulus Interval', 'type': 'float', 'value': 1., 'limits': [0.2, 300.], 
                    'suffix': 's', 'default': 1.0, 'tip': 'Time between stimuli in a sweep'},
                {'name': 'Randomize', 'type': 'bool', 'value': False, 'default': False, 'tip': 'Randomize presentation order in all dimensions'},
                {'name': 'Duration', 'type': 'float', 'value': 0.2, 'step': 0.05, 'limits': [0.001, 10],
                    'suffix': 's', 'default': 0.2, 'tip': 'Sound duration, in seconds'},
                {'name': 'Delay', 'type': 'float', 'value': 0.01, 'step': 0.05, 'limits': [0.001, 10],
                    'suffix': 's', 'default': 0.2, 'tip': 'Sound delay, in seconds'},
             ]},
             {'name': 'Clicks', 'type': 'group', 'children': [
                  {'name': 'Interval', 'type': 'float', 'value': 50., 'step': 5.0, 'limits': [1., 1000.0],
                    'suffix': 'ms', 'default': 50.0},
                  {'name': 'Number', 'type': 'int', 'value': 4, 'step': 1, 'limits': [1, 200.0],
                    'default': 4},
                  {'name': 'Duration', 'type': 'float', 'value': 1e-4, 'step': 10e-6, 
                      'limits': [10e-6, 1e-3], 'suffix': 's', 'default': 1e-4},
             ]},

             {'name': 'Modulation/CMMR', 'type': 'group', 'children': [
                  {'name': 'Frequency', 'type': 'float', 'value': 40., 'step': 5.0, 'limits': [1., 1000.0],
                    'suffix': 'Hz', 'default': 40.0},
                  {'name': 'Depth', 'type': 'float', 'value': 50.0, 'step': 5.0, 'limits': [0.0, 200.0],
                    'suffix': '%', 'default': 50.},
                  {'name': 'CMMR Flanking Type', 'type': 'list', 'values': ['None', '3Tone', 'Noise'], 'value': 'Noise'},
                  {'name': 'CMMR Flanking Phase', 'type': 'list', 'values': ['Comodulated', 'Deviant', 'Random'],
                   'value': 'Noise'},
                  {'name': 'CMMR Flanking Bands', 'type': 'int', 'value': 3, 'step': 1, 'limits': [0, 10],
                   'default': '3'},
                  {'name': 'CMMR Flanking Spacing', 'type': 'float', 'value': 0.5, 'step': 1/8., 
                      'limits': [1/16., 2.], 'suffix': 'octaves', 'default': 0.5},
             ]},
             {'name': 'RSS Params', 'type': 'group', 'children': [
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
                    

            {'name': 'File From Disk', 'type': 'str', 'value': 'test.wav', 'default': 'test.wav'},
            
            # {'name': 'Show Waveform', 'type': 'action', 'height': 12},
            # {'name': 'Show Spectrum', 'type': 'action'},
            # {'name': 'Run', 'type': 'action'},
            # {'name': 'Pause', 'type': 'action'},
            # {'name': 'Stop', 'type': 'action'},
        ]
        
        self.ptree = ParameterTree()
        self.ptreedata = Parameter.create(name='params', type='group', children=params)
        self.ptree.setParameters(self.ptreedata)
        ptreewidth = 120
        
        # now build the ui
        # hardwired buttons
        self.btn_waveform = QtGui.QPushButton("Show Waveform")
        self.btn_spectrum = QtGui.QPushButton("Show Spectrum")
        self.btn_run = QtGui.QPushButton("Run")
        self.btn_pause = QtGui.QPushButton("Pause")
        self.btn_continue = QtGui.QPushButton("Continue")
        self.btn_stop = QtGui.QPushButton("Stop")
        self.btn_quit = QtGui.QPushButton("Quit")
        
        hbox = QtGui.QGridLayout()
        hbox.setColumnStretch(0, 1)
        hbox.setColumnStretch(1, 1)
        
        hbox.addWidget(self.btn_waveform, 0, 0, 1, 1)
        hbox.addWidget(self.btn_spectrum, 0, 1, 1, 1)
        hbox.addWidget(self.btn_run, 1, 0, 1, 1)
        hbox.addWidget(self.btn_pause, 1, 1, 1, 1)
        hbox.addWidget(self.btn_continue, 1, 2, 1, 1)
        hbox.addWidget(self.btn_stop, 0, 2, 1, 1)
        hbox.addWidget(self.btn_quit, 0, 3, 1, 1)
        
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


        self.img = pg.ImageView() # view=self.plots['Spec'])
        arr = np.random.random((100, 32))
        self.img.setImage(arr)
        self.img.ui.roiBtn.hide()
#        self.img.ui.menuBtn.hide()
        self.img.show()

        glayout.nextRow()
        l2 = glayout.addLayout(colspan=3, border=(50,0,0))  # embed a new layout
        l2.setContentsMargins(10,10,10,10)
        self.plots['Plot1'] = l2.addPlot(Title="Plot1")
#        self.l2.addWidget(self.plots['Plot1'])
        self.plots['Plot1'].getAxis('bottom').setLabel('F (kHz)')
        self.plots['Plot1'].getAxis('left').setLabel('dB SPL')
        self.plots['Plot1'].setTitle('FRA')
        self.plots['Plot1'].setXRange(2, 72, padding=0)
        #self.plots['Plot1'].setLogMode(x=True)
        self.plots['Plot1'].setYRange(-10, 110, padding=0)
        xd = np.arange(4, 64, 1)
        yd = np.arange(0, 100, 2.5)
        spots = []
        self.lastPoint = None
        for i in range(xd.shape[0]):
            for j in range(yd.shape[0]):
                spots.append({'pos': (xd[i], yd[j]), 'size': 5, 'pen': {'color': 'w', 'width': 1, 'alpha': 0.5},
                    'brush': pg.mkBrush('b')})
        self.spi = pg.ScatterPlotItem(size=5, pen=pg.mkPen('k'), brush=pg.mkBrush('b'), symbol='s')
        self.spi.addPoints(spots)
        self.plots['Plot1'].addItem(self.spi)
        self.spi.sigClicked.connect(self.getClickedLocation)
        
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
        self.controller = Controller(self.ptreedata, self.plots, self.img)

        self.controller.setAllParameters(params)  # synchronize parameters with the tree
        self.ptreedata.sigTreeStateChanged.connect(self.controller.change)  # connect parameters to their updates
        
        # now connect the buttons
        self.btn_waveform.clicked.connect(self.controller.show_wave)
        self.btn_spectrum.clicked.connect(self.controller.show_spectrogram)
        self.btn_run.clicked.connect(self.controller.run)
        self.btn_pause.clicked.connect(self.controller.pause_run)
        self.btn_continue.clicked.connect(self.controller.continue_run)
        self.btn_stop.clicked.connect(self.controller.stop_run)
        self.btn_quit.clicked.connect(self.controller.quit)

    def getClickedLocation(self, points):
        # print (dir(points))
        # print (points.event())
        # print('mouse click: ', points.mouseClickEvent(points.event()))
        # print('mouse doubleclick: ', points.mouseDoubleClickEvent(points.event()))
        self.mousePoint = points.ptsClicked[0].viewPos()
        print('mousepoint: ', self.mousePoint.x(), self.mousePoint.y()) 
        points.ptsClicked[0].setBrush(pg.mkBrush('r'))
        points.ptsClicked[0].setSize(5)
        if self.lastPoint is None:
            self.lastPoint = points.ptsClicked[0]
        else:
            self.lastPoint.setBrush(pg.mkBrush('b'))
            self.lastPoint.setSize(5)
            self.lastPoint = points.ptsClicked[0]
        self.controller.protocol = 'Tone RI'
        self.controller.tone_frequency = self.mousePoint.x()
        self.controller.dbspl = self.mousePoint.y()
        self.controller.prepare_run()
        self.controller.run()
    
if __name__ == '__main__':
    gui = BuildGui()
    
    ## Start Qt event loop unless running in interactive mode.
    ## Event loop will wait for the GUI to activate the updater and start sampling.
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
    gui.controller.quit()