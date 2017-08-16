"""
Read frequency.cal file from matlab abr program calibration, transate results to dictionary
"""
import scipy.io
import os.path
import numpy as np
import scipy.interpolate
from collections import OrderedDict
import matplotlib.pyplot as mpl


class CalReader():
	def __init__(self, filename=None):
		if filename is None:
			filename = 'C:/Users/experimenters/Desktop/ABR_Code/frequency_MF1.cal'
		print('Calreader from : %s', filename)
		d = scipy.io.loadmat(filename)

		self.Freqs = d['CAL']['Freqs'].item().ravel()  # frequencies tested
		self.maxdb = d['CAL']['maxdB'].item().ravel()  # db measured re max
		self.refspl = d['CAL']['RefSPL'].item().ravel()[0]  # reference for max (typically, 110)
		self.calattn = d['CAL']['CalAttn'].item().ravel()[0]
		self.caldate = d['CAL']['Date'].item().ravel()[0]
		self.finterp = scipy.interpolate.interp1d(self.Freqs, self.maxdb, kind='cubic')


	def plotcal(self):

		ax = mpl.subplot(111)
		ax.semilogx(self.Freqs, self.maxdb, 'ro-')

	def getSPLatF(self, f, spl):
		#splatf = np.interp(f, self.Freqs, self.maxdb)
		splatf = self.finterp(f)
		attn = splatf + self.calattn - spl

		print('F: %f   splatf: %f   desired spl: %f   attn: %f' % (f, splatf, spl, attn))
		return attn

if __name__ == "__main__":
	c = CalReader()
	c.plotcal()
	c.getSPLatF(16000, 75.1)
	mpl.show()
