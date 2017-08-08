from __future__ import print_function

import sys
import datetime
import pickle
import platform
import serial
import time
import numpy as np
import scipy.signal
import pyqtgraph as pg
from PyQt4 import QtGui, QtCore
from pyqtgraph.parametertree import Parameter, ParameterTree
import pystim

class Controller(object):
    def __init__(self):
        self.toneFrequency = 4000.  # tone pip frequency
        self.RI = '20;100/10'  # seqparse for Rate-Intensity intensity series
        self.FRA = '4;48/8l'  # seqparse for tone pip FRA
        self.dMod = 0.  # modulation depth
        self.fMod = 40. # modulation frequency
        self.RF = 2.5e-3  # rise-fall time
        
        self.randomize = False  # randomize order of presentation (or not...)
        
        
        pass
        
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
        pass
    

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
            path = ptreedata.childPath(param)
            # if path is not None:
            #     childName = '.'.join(path)
            # else:
            #     childName = param.name()
            #
            # Parameters and user-supplied information
            #
            if path[1] == 'tone':
                self.filename = data
            if path[1] == 'Interval':
                self.readInterval = data
            if path[1] == 'Duration':
                self.readDuration = data
            if path[1] == 'Invert':
                self.invertData = data
            if path[1] == 'MaxSamples':
                self.maxSamples = data
            if path[1] == 'LPF':
                self.LPFFreq = data
            if path[1] == 'Notch':
                self.NotchFreq = data
            if path[1] == 'NotchEnabled':
                self.NotchEnabled = data
            if path[1] == 'Info':
                self.InfoText = data
            #
            # Actions:
            #
            if path[1] == 'Start New':
                self.startRun()
            if path[1] == 'Stop/Pause':
                self.stopRun()
            if path[1] == 'Continue':
                self.continueRun()
            if path[1] == 'Save Visible':
                self.storeData()
            if path[1] == 'Load File':
                fn = self.getFilename()
                if fn is not None:
                    self.loadData(filename=fn)
            if path[1] == 'New Filename':
                self.filename = self.makeFilename()
    
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
            if p == 'Filename':
                self.filename = p['value']
            if p == 'Interval':
                self.readInterval = p['value']
            if p == 'Duration':
                self.readDuration = p['value']
            if p == 'Invert':
                self.invertData = p['value']
            if p == 'MaxSamples':
                self.maxSamples = p['value']
            if p == 'LPF':
                self.LPFFreq = p['value']
            if p == 'Notch':
                self.NotchFreq = p['value']
            if p == 'NotchEnabled':
                self.NotchEnabled = p['value']
            if p == 'Info':
                self.InfoText = p['value']

    def startRun(self):
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
        
        self.runtime = 0
        self.NSamples = 0
        self.startTime = datetime.datetime.now()
        self.prepareRun()  # reset the data arrays
        self.continueRun()
    
    def continueRun(self):
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
        self.timer.timeout.connect(updater.update)
        self.timedWrite = QtCore.QTimer()
        self.timedWrite.timeout.connect(updater.storeData)
        self.update() # do the first update, then start time
        self.timer.start(self.readInterval * 1000)
        self.timedWrite.start(30 * 1000.)  # update file in 1 minute increments
    
    def stopRun(self):
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
        self.timedWrite.stop()
        self.storeData()

    def prepareRun(self):
        """
        Clear out all arrays for the data collection run
        
        Parameters
        ----------
        None
        
        Returns
        -------
        Nothing
        """
        self.runningRate = []
        self.runningVar = []
        self.runningTime = []
        self.RRInterval = []
        self.varplot = None
        self.currentWave = []
        self.out = []
        self.makeFilename()
        
    def setLPF(self, freq):
        """
        Store the low-pass filter frequency to be used in analysis
        
        Parameters
        ----------
        freq: float, cutoff frequency in Hz (no default)
        
        Returns
        -------
        Nothing
        """
        self.ecg.setLPF(freq)
        self.LPFFreq = freq
        
    def setNotch(self, freq):
        """
        Store the Notch filter frequency to be used in analysis
    
        Parameters
        ----------
        freq: float, Notch frequency in Hz (no default)
    
        Returns
        -------
        Nothing
        """
        self.ecg.setNotch(freq)
        self.NotchFreq = freq
# Build GUI and window

class BuildGui(object):
    def __init__(self):
        app = pg.mkQApp()
        self.win = QtGui.QWidget()
        self.layout = QtGui.QGridLayout()
        self.win.setLayout(self.layout)
        self.win.show()
        self.win.setWindowTitle('Stim Controller')
        self.win.setGeometry( 100 , 100 , 1024 , 600)
        #win.resize(1024,800)

        # Define parameters that control aquisition and buttons...
        params = [
            {'name': 'Stimulus', 'type': 'group', 'children': [
                {'name': 'Repetitions', 'type': 'int', 'value': 10, 'limits': [1, 10000], 'default': 10},
                {'name': 'Interstimulus Interval', 'type': 'float', 'value': 1., 'limits': [0.2, 300.], 
                    'suffix': 's', 'default': 1.0},
                {'name': 'Duration', 'type': 'float', 'value': 1., 'step': 0.5, 'limits': [0.5, 10],
                    'suffix': 's', 'default': 0.2},
                {'name': 'Tone Frequency', 'type': 'float', 'value': 4.0., 'step': 1.0, 'limits': [0.5, 99.0],
                    'suffix': 'kHz', 'default': 4.0},
                {'name': 'Rise-Fall', 'type': 'float', 'value': 2.5, 'step': 0.5, 'limits': [0.5, 20.],
                    'suffix': 'ms', 'default': 2.5},
                {'name': 'Intensities', 'type': 'str', 'value': '20;90/10',
                    'suffix': 'dB', 'default': '20;90/10'},
                {'name': 'Duration', 'type': 'float', 'value': 1., 'step': 0.5, 'limits': [0.5, 10],
                    'suffix': 's', 'default': 0.2},
                    

                {'name': 'Protocol', 'type': 'list', 'values': ['Tone RI', 'Noise RI', 'FRA',
                        'RSS', 'DMR', 'SSN', 'Tone SAM', 'Noise SAM', 'Clicks',
                        'NotchNoise', 'BandPass Noise'], 'value': 'Tone RI'},

                {'name': 'File From Disk', 'type': 'str', 'value': 'test.wav', 'default': 'test.wav'},
            
                {'name': 'Show Waveform', 'type': 'action', 'height': 12},
                {'name': 'Show Spectrum', 'type': 'action'},
                {'name': 'Run', 'type': 'action'},
                {'name': 'Pause', 'type': 'action'},
                {'name': 'Stop', 'type': 'action'},
        ]}]
        ptree = ParameterTree()
        ptreedata = Parameter.create(name='params', type='group', children=params)
        ptree.setParameters(ptreedata)
        print (dir(ptree))
        exit(1)# for child in ptree.children():
        #     print ('child: ', child)
        #  #   print (dir(child)
        #     if isinstance(child, QtGui.QWidget):
        #         child.setMaximumHeight(20)
        #     for ch2 in child.children():
        #         print('child2: ', ch2)
        #         for ch3 in ch2.children():
        #             print ('ch3: ', ch3)
        #             print (dir(ch3))
        #             for ch4 in ch3.children():
        #                 print ('ch4: ', ch4)
        #                 print (dir(ch4))
        #     #     print (dir(ch2))
            #     ch2.setMaximumHeight(12)
        # build layout for plots and parameters
        self.layout.addWidget(ptree, 0, 0, 5, 2) # Parameter Tree on left

        # add space for the graphs
        view = pg.GraphicsView()
        l = pg.GraphicsLayout(border=(50,50,50))
        view.setCentralItem(l)
        self.layout.addWidget(view, 0, 3, 5, 3)  # data plots on right
        self.plots = {}
        self.plots['Wave'] = l.addPlot()
        self.plots['Wave'].getAxis('left').setLabel('V', color="#ff0000")
        self.plots['Wave'].setTitle('Waveform', color="#ff0000")
        self.plots['Wave'].getAxis('bottom').setLabel('t (min)', color="#ff0000")
        self.plots['Wave'].setYRange(0, 300)
        ## create a new ViewBox, link the right axis to its coordinate system
        l.nextRow()
        self.plots['Spec'] = l.addPlot() # 
        self.plots['Spec'].setXLink(self.plots['Wave'])
        self.plots['Spec'].getAxis('left').setLabel('Intensity', color='#0000ff')
        self.plots['Spec'].setTitle('Freq', color="#ff0000")
        self.plots['Spec'].setYRange(0, 20)
        self.plots['Spec'].getAxis('bottom').setLabel('t (min)', color='#0000ff')

        l.nextRow()
        l2 = l.addLayout(colspan=3, border=(50,0,0))  # embed a new layout
        l2.setContentsMargins(10,10,10,10)
        self.plots['Plot1'] = l2.addPlot(Title="Plot1")
        self.plots['Plot1'].getAxis('bottom').setLabel('t (s)')
        self.plots['Plot1'].getAxis('left').setLabel('V')
        self.plots['Plot1'].setTitle('Plot 1')
        self.plots['Plot2'] = l2.addPlot(Title="Plot2")
        self.plots['Plot2'].getAxis('bottom').setLabel('t (s)')
        self.plots['Plot2'].getAxis('left').setLabel('V')
        self.plots['Plot2'].setTitle('Plot 2')
        self.plots['Plot3'] = l2.addPlot(Title="Plot 3")
        self.plots['Plot3'].setTitle('Plot3')
        self.plots['Plot3'].getAxis('bottom').setLabel('t (s)')
        self.plots['Plot3'].getAxis('left').setLabel('V')


        #
        # Initialize the updater with needed information about the plots
        #
        controller = Controller()

        controller.setAllParameters(params)  # synchronize parameters with the tree

        ptreedata.sigTreeStateChanged.connect(controller.change)  # connect parameters to their updates

if __name__ == '__main__':
    gui = BuildGui()
    
    ## Start Qt event loop unless running in interactive mode.
    ## Event loop will wait for the GUI to activate the updater and start sampling.
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()