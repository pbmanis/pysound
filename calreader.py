"""
Read frequency.cal file from matlab abr program calibration, transate results to dictionary
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
		if filename is None:
			filename = 'C:/Users/experimenters/Desktop/ABR_Code/frequency_MF1.cal'
		self.cals = OrderedDict()
		self.plots = None
		#self.readfile(filename)

	def readfiles(self):
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
		print('Calreader from : %s', filename)
		d = scipy.io.loadmat(filename)

		self.cals[filename] = OrderedDict()
		self.cals[filename]['frequencies'] = d['CAL']['Freqs'].item().ravel()  # frequencies tested
		self.cals[filename]['maxdb'] = d['CAL']['maxdB'].item().ravel()  # db measured re max
		self.cals[filename]['refspl'] = d['CAL']['RefSPL'].item().ravel()[0]  # reference for max (typically, 110)
		try:  # missing field in early versions
			self.cals[filename]['calattn'] = d['CAL']['CalAttn'].item().ravel()[0]
		except:
			self.cals[filename]['calattn'] = 30.
		self.cals[filename]['caldate'] = d['CAL']['Date'].item().ravel()[0]
		self.cals[filename]['finterp'] = scipy.interpolate.interp1d(self.cals[filename]['frequencies'],
										self.cals[filename]['maxdb'], kind='cubic')

	def plotcals(self):
		for filename in self.cals.keys():
		#	print self.cals[filename]
			if self.plots is None:
				self.plots = mpl.subplot(111)
			h, shortfn = os.path.split(filename)
			self.plots.semilogx(self.cals[filename]['frequencies'], self.cals[filename]['maxdb'], label=shortfn)
		self.plots.grid(True, which='both')
		mpl.legend(loc=3,
           ncol=2, mode="expand", borderaxespad=0.)

	def getSPLatF(self, f, spl, filename):
		#splatf = np.interp(f, self.Freqs, self.maxdb)
		splatf = self.cals[filename]['finterp'](f)
		attn = splatf + self.cals[filename]['calattn'] - spl

		print('F: %f   splatf: %f   desired spl: %f   attn: %f' % (f, splatf, spl, attn))
		return attn

if __name__ == "__main__":
	app = QtGui.QApplication(sys.argv)
	c = CalReader()
	c.readfiles()
	c.plotcals()
	c.getSPLatF(16000, 75.1, c.cals.keys()[0])
	mpl.show()
