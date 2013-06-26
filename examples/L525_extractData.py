import ixppy
import pylab as p
import numpy as np
import copy
import sys
h5 = ixppy.toolsHdf5
h5w = ixppy.toolsHdf5.datasetWrite
ttfit = ixppy.toolsTiming.findStepWithErfFit 

roiTT    = {'limits': np.array([ 146.,  198.]), 'projection': 'vertical range'}
roiTTref = {'limits': np.array([ 63. ,  105.]), 'projection': 'vertical range'}

def extractData(run,calibcycles=None,Nmax=None):
	f=open("run%03d_extract.txt" % run,"w")
	sys.stderr = f
	sys.stout  = f
	dets = ['ipm2','ipm3',"diodeU","diode3","timeTool","opal0"]
	d=ixppy.dataset(('xpp52512',run),dets)
	d.filtTimestamps()
	ixppy.TTextractProfiles(d.opal0,d.ipm2.sum,Nmax=Nmax,refthreshold=0.01,profileLimits=roiTT,profileLimitsRef=roiTTref,calibcycles=calibcycles)

	h=h5.openOrCreateFile("run%03d_extract.h5" % run,"w")
	h.attrs["ScanMot"] = d.scanMot
	h.attrs["ScanVec"] = d.scanVec
	my_steps = []
	x=np.arange(1024)
	n=200
	N=800
	fitrange=100
	h.attrs["TTfitinfo"] = "from %s to %s (fitrange %d)" % (n,N,fitrange)
	for i in range(len(d.opal0.TTtraces)):
		my_stepsc=[]
		for s in range(len(d.opal0.TTtraces[i])):
			if (d.ipm2.sum[i][s] < 0.01):
				my_stepsc.append([0,0,0])
			else:
				y = d.opal0.TTtraces[i][s][n:N]
				(step,amp,sig,xfit,yfit)=ttfit(x[n:N],y,fitrange=fitrange)
				my_stepsc.append( [step,amp,sig] )
		my_steps.append ( np.array(my_stepsc) )
	for i in range(len(d.opal0.TTtraces)):
		h5w(h,"/c%03d/TTtraces" % i,d.opal0.TTtraces[i].value)
		h5w(h,"/c%03d/ipm2" % i,d.ipm2.data[i])
		h5w(h,"/c%03d/ipm3" % i,d.ipm3.data[i])
		h5w(h,"/c%03d/diodeU" % i,d.diodeU.data[i])
		h5w(h,"/c%03d/diode3" % i,d.diode3.data[i])
		h5w(h,"/c%03d/TTXPP/pos" % i,d.timeTool.fltpos[i])
		h5w(h,"/c%03d/TTXPP/ampl" % i,d.timeTool.ampl[i])
		h5w(h,"/c%03d/TTXPP/fwhm" % i,d.timeTool.fltposfwhm[i])
		h5w(h,"/c%03d/TTERF/pos" % i,my_steps[i][:,0])
		h5w(h,"/c%03d/TTERF/amp" % i,my_steps[i][:,1])
		h5w(h,"/c%03d/TTERF/fwhm" % i,my_steps[i][:,2]*2.35)
	f.close()
	sys.exit(0)

if (__name__=="__main__"):
	import sys
	run = int(sys.argv[1])
	ixppy.tools.removeFileIfExists("/home/marcoc/ixppy_cache/xpp52512-r%04d.h5" % run)
	ixppy.tools.removeFileIfExists("/reg/neh/home/marcoc/ixppy_cache/xpp52512-r%04d.h5" % run)
	ixppy.tools.removeFileIfExists("run%03d_extract.h5" % run)
	extractData(run)
