from future import __print_function__
"""
Read frequency.cal file(s) from matlab abr program calibration.

Provides:
    Reading of the file and storage in an instance of the CalReader Class.
    Plotting of calibration data (including multiple files for comparison).
    Interpolation of calibration curve to return attenuation required
    to set the SPL at a single frequency.

TODO:
    Provide spectral correction of complex (non-tone) waveform - flattening,
    based on calibration.

8/16/2017 pbm

"""

import sys
import scipy.io
import os.path
import numpy as np
import scipy.interpolate
from collections import OrderedDict
import matplotlib.pyplot as mpl
from PyQt4 import QtGui, QtCore


class CalReader():
    def __init__(self, filename=None):
        """
        CalReader reads matlab calibration files created by the 'ABR4' program
        
        Parameters
        ----------
        filename : string (default: None)
            name of file to read; if None, then just set up
        """
        self.cals = OrderedDict()  # storage of all calibrations read
        self.plots = None
        
        if filename is None:
            self.readfiles() # use a gui to grab a calibration file (or more)
#        filename = 'C:/Users/experimenters/Desktop/ABR_Code/frequency_MF1.cal'
        #self.readfile(filename)
    
    def readfiles(self):
        """
        Read calibration files into the calibration dictionary. Opens a Qt file dialo
        
        Parameters
        ----------
        None
        """
        
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.ExistingFiles)
        dlg.setFilter("Matlab ABR Calibration files (*.cal)")
        #self.filenames = dlg(None, "Open files", "C:\Users\Desktop\ABR_Code", "Cal files (*.cal)")
        #print [str(fn) for fn in self.filenames]
        
        if dlg.exec_():
            self.filenames = dlg.selectedFiles()
            print('filenames: ', self.filenames)
            for f in self.filenames:
                print str(f)
                self.readfile(str(f))
    
    def readfile(self, filename):
        """
        Read one file and store the results in the cals dictionary using
        the filename as the key
        
        Parameters
        ----------
        None
        """

#        print('Calreader from : %s', filename)
        d = scipy.io.loadmat(filename)
        
        self.cals[filename] = OrderedDict()
        self.cals[filename]['frequencies'] = d['CAL']['Freqs'].item().ravel()  # frequencies tested
        self.cals[filename]['maxdb'] = d['CAL']['maxdB'].item().ravel()  # db measured re max
        self.cals[filename]['refspl'] = d['CAL']['RefSPL'].item().ravel()[0]  # reference for max (typically, 110)
        try:  # handle the missing calattn field in early versions
            self.cals[filename]['calattn'] = d['CAL']['CalAttn'].item().ravel()[0]
        except:
            self.cals[filename]['calattn'] = 30.
        self.cals[filename]['caldate'] = d['CAL']['Date'].item().ravel()[0]
        self.cals[filename]['finterp'] = scipy.interpolate.interp1d(self.cals[filename]['frequencies'],
                                        self.cals[filename]['maxdb'], kind='cubic')
    
    def plotcals(self):
        """
        Plot the calibration data for multiple files
        
        Parameters
        ----------
        None
        """
        
        for filename in self.cals.keys():
        #    print self.cals[filename]
            if self.plots is None:
                self.plots = mpl.subplot(111)
            h, shortfn = os.path.split(filename)
            self.plots.semilogx(self.cals[filename]['frequencies'], self.cals[filename]['maxdb'], label=shortfn)
        self.plots.grid(True, which='both')
        mpl.legend(loc=3,
           ncol=2, mode="expand", borderaxespad=0.)
    
    def getAttnforSPL(self, f, spl, filename):
        """
        Return the attenuator setting needed to get the requested SPL using
        the calibration from filename
        
        Parameters
        ----------
        f : float, Hz (no default)
            Frequency to request the attenuation
        spl : float, dB SPL (no default)
            desired sound pressure level
        filename : str (no default)
            calibration file to use for interpolation
        
        Returns
        -------
        float : attn (dB)
        
        """
        
        #splatf = np.interp(f, self.Freqs, self.maxdb)
        splatf = self.cals[filename]['finterp'](f)
        attn = splatf + self.cals[filename]['calattn'] - spl
        
#        print('F: %f   splatf: %f   desired spl: %f   attn: %f' % (f, splatf, spl, attn))
        return attn

if __name__ == "__main__":
    """
    Testing and plotting multiple calibrations
    """
    app = QtGui.QApplication(sys.argv)
    c = CalReader()
    c.readfiles()
    c.plotcals()
    c.getAttnforSPL(16000, 75.1, c.cals.keys()[0])  # just to check result matches checkcal
    mpl.show()
