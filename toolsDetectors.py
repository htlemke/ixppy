import numpy as np
import pylab as pl
import toolsPlot

def addScanVecToSingleShotReadings(scanv,tt):
	""" tt must be either a list of vectors or a matrix """
	if isinstance(tt,list):
		return [scanv[i]+tt[i][:] for i in range(len(scanv))]
	elif (tt.shape[0] == len(scanv)):
		return [scanv[i]+tt[i,:] for i in range(len(scanv))]
	elif (tt.shape[1] == len(scanv)):
		return [scanv[i]+tt[:,i] for i in range(len(scanv))]

class nonLinearCorrection(object):
	""" Class to be used for correcting detector non linearity
	The reason why is a class and not a simple fucntion is to be able
	to apply calculate the correction factors for certain calibcycles
	(like without laser or for negative time delays) and then use the
	pre-calculated factors to correct for any time delay
	usage:
	nonlin = nonLinearCorrection()
	nonlin.calibrate(det[calib1],mon[calib1])
	for c in calis:
		(mon[c],det[c])=nonlin.correct(mon[c],det[c])
		.....
	"""
	def __init__(self):
		pass

	def calibrate(self,mon,dio,plot=False):
		self.mon = mon
		self.dio = dio
		self.poly = np.polyfit(mon,dio,2)
		if (plot):
			toolsPlot.nfigure("calibration")
			ax1=pl.subplot("211",title="det vs monitor (before correction)")
			pl.plot(mon,dio,"+")
			m=np.min(mon); M=np.max(mon)
			x = np.arange(m,M,(M-m)/100.)
			pl.plot(x,np.polyval(self.poly,x))
			pl.subplot("212",title="Ratio vs monitor",sharex=ax1)
			#pl.plot(mon,dio/mon,"+",label="before")
			M=mon.max()
			m=mon.min()
			h1 = np.histogram(mon,np.arange(m,M,(M-m)/30.),weights=dio/mon)
			h0 = np.histogram(mon,np.arange(m,M,(M-m)/30.),)
			pl.plot(toolsPlot.histVecCenter(h1[1]),h1[0]/h0[0],label="before")
			mon,dio=self.correct(mon,dio)
			h1 = np.histogram(mon,np.arange(m,M,(M-m)/30.),weights=dio/mon)
			h0 = np.histogram(mon,np.arange(m,M,(M-m)/30.),)
			pl.plot(toolsPlot.histVecCenter(h1[1]),h1[0]/h0[0],label="after")
			pl.legend()

	def correct(self,mon,dio,plot=False):
		pp = self.poly
		x0 = (-pp[1]+np.sqrt(pp[1]**2-4*pp[0]*pp[2]))/2/pp[0]
		if ( x0<0 ):
			dioC	= dio-pp[2]
			monC	= mon
		else:
			monC	= mon-x0
			dioC	= dio
		if (pp[0]>0):
			monC	= monC + pp[0]/pp[1]*monC**2
		else:
			dioC	= dioC - pp[0]*monC**2
		if (plot):
			toolsPlot.nfigure("Check non lin correction")
			pl.plot(mon,dio/mon,"+",label = "before correction")
			pl.plot(monC,dioC/monC,"o",label = "after correction")
			pl.ylabel("ratio signal/monitor")
			pl.legend()
			pl.grid()
		return (monC,dioC)

