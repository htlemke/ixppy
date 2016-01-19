from __future__ import print_function
import numpy as np
import pylab as plt
import toolsPlot

def corrNonlin(data,polypar,data_0=0,correct_0=0):
  """ unses parameters found by corrNonlinGetPar to correct data;
  example of usage (assuming mon is non linear and loff is
  the laseroff filter
  d = ixppy.dataset("xppc3614-r0102.stripped.h5")
  loff = d.eventCode.code_91.filter(True)
  mon = d.ipm3.sum;
  dio = d.diodeU.channel0;
  poly = ixppy.tools.corrNonlinGetPar(dio*loff,mon*loff,plot=True)
  mon_corr = ixppy.corrNonlin(mon,poly)"""
  m = 1/np.polyval(np.polyder(polypar),data_0)
  return m*(np.polyval(polypar,data)-correct_0) + data_0
  
def corrNonlinGetPar(data,correct,order=2,data_0=0,correct_0=0,
    displayWarning=True,plot=False):
  """ Find parameters for non linear correction
    *data* should be an 1D array (use .ravel() in case) of the
    detectors that is suspected to be non linear
    *correct* is the detector that is sussposed to be linear
    *data_0" is an offset to use for the data (used only if plotting"
    *correct_0* offset of the "linear detector"""
  # poor man wrapping :D #
  try:
    data = data.ravel()
  except AttributeError:
    pass
  try:
    correct = correct.ravel()
  except AttributeError:
    pass
  p =  np.polyfit(data,correct,order)
  if order>=2 and p[-3]<0:
    print("corrNonlinGetPar: consistency problem, second order coefficient should \
    be > 0, please double check result (plot=True) or try inverting the data and the\
    correct arguments")
  p[-1] = p[-1]-correct_0
  if plot:
    d = corrNonlin(data,p,data_0=data_0,correct_0=correct_0)
    plt.plot(correct,data,".",label="before correction")
    plt.plot(correct,d,".",label="after correction")
    poly_lin = np.polyfit(correct,d,1)
    xmin = min(correct.min(),0)
    xtemp = np.asarray( (xmin,correct.max()) )
    plt.plot(xtemp,np.polyval(poly_lin,xtemp),
       label="linear fit")
    plt.plot(correct,d-np.polyval(poly_lin,correct),
       ".",label="difference after-linear")
    plt.xlabel("correct")
    plt.ylabel("data")
    plt.legend()
  return p

def addScanVecToSingleShotReadings(scanv,tt):
  """ tt must be either a list of vectors or a matrix; now it
  is not needed as mem data natively supports that:
  d.timeTool.pos + d.scan.lxt"""
  print("this funciton is obsolete: please use d.timeTool.pos + d.scan.lxt")
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
    print("This funciton is obsolete, please use: corrNonlinGetPar and corrNonlin")
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
      dioC  = dio-pp[2]
      monC  = mon
    else:
      monC  = mon-x0
      dioC  = dio
    if (pp[0]>0):
      monC  = monC + pp[0]/pp[1]*monC**2
    else:
      dioC  = dioC - pp[0]*monC**2
    if (plot):
      toolsPlot.nfigure("Check non lin correction")
      pl.plot(mon,dio/mon,"+",label = "before correction")
      pl.plot(monC,dioC/monC,"o",label = "after correction")
      pl.ylabel("ratio signal/monitor")
      pl.legend()
      pl.grid()
    return (monC,dioC)

