from __future__ import print_function

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

from pathlib import Path
import sys
from collections import OrderedDict
from typing import Union
from tkinter import Tk
import tkinter.filedialog as FD

import numpy as np
import scipy.interpolate
import scipy.io
import matplotlib
matplotlib.use('Qt5Agg')
from matplotlib import pyplot as mpl


class CalReader:
    def __init__(self, filename:Union[str, None]=None):
        """
        CalReader reads matlab calibration files created by the 'ABR4' program
        
        Parameters
        ----------
        filename : string (default: None)
            name of file to read; if None, then just set up
        """
        # self.app = QtGui.QApplication(sys.argv)
        
        self.cals = OrderedDict()  # storage of all calibrations read
        self.plots = None
        self.filenames = None  # or list
        if filename is None:
            self.readfiles()  # use a gui to grab a calibration file (or more)

    #        filename = 'C:/Users/experimenters/Desktop/ABR_Code/frequency_MF1.cal'
    # self.readfile(filename)

    def readfiles(self):
        """
        Read calibration files into the calibration dictionary. Opens a Qt file dialo
        
        Parameters
        ----------
        None
        """
       
        window = Tk()
        window.wm_attributes('-topmost', 1)
        window.withdraw()   # this supress the tk window
        filenames = FD.askopenfilenames(parent=window, title="Get Calibration Files", 
            initialdir="../tests/test_cal_data", filetypes=(("Cal", "*.cal"),),
            multiple=True)
        self.filenames = filenames
        print(self.filenames)
        for f in self.filenames:
            print("File: {:s}".format(str(f)))
            self.readfile(str(f))
        window.destroy()


    def readfile(self, filename:str):
        """
        Read one file and store the results in the cals dictionary using
        the filename as the key
        
        Parameters
        ----------
        filename : the name of the calibration file to read.
        
        """

        #        print('Calreader from : %s', filename)
        d = scipy.io.loadmat(filename)

        self.cals[filename] = OrderedDict()
        self.cals[filename]["frequencies"] = (
            d["CAL"]["Freqs"].item().ravel()
        )  # frequencies tested
        self.cals[filename]["maxdb"] = (
            d["CAL"]["maxdB"].item().ravel()
        )  # db measured re max
        self.cals[filename]["refspl"] = (
            d["CAL"]["RefSPL"].item().ravel()[0]
        )  # reference for max (typically, 110)
        try:  # handle the missing calattn field in early versions
            self.cals[filename]["calattn"] = d["CAL"]["CalAttn"].item().ravel()[0]
        except:
            self.cals[filename]["calattn"] = 30.0
        self.cals[filename]["caldate"] = d["CAL"]["Date"].item().ravel()[0]
        self.cals[filename]["finterp"] = scipy.interpolate.interp1d(
            self.cals[filename]["frequencies"],
            self.cals[filename]["maxdb"],
            kind="cubic",
        )

    def plotcals(self):
        """
        Plot the calibration data for multiple files
        
        Parameters
        ----------
        None
        """

        for filename in list(self.cals.keys()):
            if self.plots is None:
                self.plots = mpl.subplot(111)
            shortfn = Path(filename).name
            self.plots.semilogx(
                self.cals[filename]["frequencies"],
                self.cals[filename]["maxdb"],
                label=shortfn,
            )
        self.plots.grid(True, which="both")
        mpl.legend(loc=3, ncol=2, mode="expand", borderaxespad=0.0)

    def getAttnforSPL(self, f:float, spl:float, filename:str):
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

        # splatf = np.interp(f, self.Freqs, self.maxdb)
        splatf = self.cals[filename]["finterp"](f)
        attn = splatf + self.cals[filename]["calattn"] - spl

        #        print('F: %f   splatf: %f   desired spl: %f   attn: %f' % (f, splatf, spl, attn))
        return attn


if __name__ == "__main__":
    """
    Testing and plotting multiple calibrations
    """
    c = CalReader()  # calls to read files

    if c.filenames is not None:
        c.plotcals()
        attn = c.getAttnforSPL(
            16000, 75.1, list(c.cals.keys())[0]
        )  # just to check result matches checkcal

        mpl.show()

