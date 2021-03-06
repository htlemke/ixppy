from __future__ import print_function
from toolsLog import logbook
import numpy as np
import scipy as sp
import numpy.linalg as linalg
from scipy.special import erf,gamma
import matplotlib
import pylab as pl
import time
import types
import numpy.ma as ma
#from toolsVarious import iterfy
#from tools import *
import ixppy

def iterfy(iterable):
    if isinstance(iterable, basestring):
        iterable = [iterable]
    try:
        iter(iterable)
    except TypeError:
        iterable = [iterable]
    return iterable


#def poiss_prob(x,count):
	#x = np.asarray(x)
	#P = np.zeros(x.shape)
	#i=0
	#for xx in x:
		#P[i] = count**xx*np*exp(-count)/sp.factorial(xx)
		#i=i+1
	#return P

def poiss_prob(x,count):
  x = np.asarray(x)
  count = float(count)
  return count**x*np.exp(-count)/gamma(x+1)

#def gauss_norm(X,xdat):
#	ydat = 1./sqrt(2.*pi*X[2]**2)*X[0]*exp(-(xdat-X[1])**2/2/X[2]**2)
#	return ydat
#
def gauss(dat,A,pos,sig):
  res = A*np.exp(-(dat-pos)**2/sig**2/2)
  return res

def photonPeaks(x,amps=None,sigma=1,gain=5):
  res = np.zeros_like(x)
  if np.iterable(amps[0]):
    ns,amps = zip(*amps)
  else:
    ns = xrange(len(amps))
  if np.iterable(sigma):
    for n,amp,sig in zip(ns,amps,sigma):
      res += gauss(x,amp,n*gain,sig)
  else:
    for n,amp in zip(ns,amps):
      res += gauss(x,amp,n*gain,sigma)
  return res

def photonPeaksPoisson(x,count=1.5,gain=5,sigma=1,fracCutoff=.001):
  assert count>=0, 'poisson Count needs to be greater or equal than zero.'
  mx = np.max(poiss_prob(np.floor(count),count),poiss_prob(np.floor(count),count))
  lolim = np.floor(count)
  while lolim >0 and poiss_prob(lolim,count)>fracCutoff*mx:
    lolim -=1
  print(lolim)
  hilim = np.ceil(count)
  while poiss_prob(hilim,count)>fracCutoff*mx:
    hilim +=1
  print(hilim)
  return photonPeaks(x,amps=zip(np.arange(lolim,hilim+1),poiss_prob(np.arange(lolim,hilim+1),count)),gain=gain,sigma=sigma)


  
  


def digitize2D(x1,x2,bins1,bins2):
	bn1 = np.digitize(x1,bins1)
	bn2 = np.digitize(x2,bins2)
	sz1 = bins1.shape[0]+1
	sz2 = bins2.shape[0]+1
	bnmat = np.reshape(arange(sz1*sz2),[sz1,sz2])
	bn = bnmat[bn1,bn2]
	return bn


def digitizeND(positions, binnings, maskInput=True):
  mask = np.zeros(np.prod(positions[0].shape),dtype=bool)
  shape = []
  inds_unravelled = []
  for pos,bins in zip(positions,binnings):
    inds = np.digitize(pos.ravel(),bins)
    if maskInput:
      #de=bug
      mask |= (inds==0)
      mask |= (inds==len(bins))
      inds -=1
      shape.append(len(bins)-1)
    else:
      shape.append(len(bins)+1)
    inds_unravelled.append(inds)
  if maskInput:
    inds_unravelled = [ti[~mask] for ti in inds_unravelled]
  inds_ravelled = np.ravel_multi_index(inds,shape)
  return inds_ravelled,shape,mask 

def bincountND(inds_ravelled,shape,mask=None,weights=None):
  if mask is not None:
    return np.bincount(inds_ravelled,
	weights=weights[mask],
	minlength=np.prod(shape)).reshape(shape)
  else:
    return np.bincount(inds_ravelled,
	weights=weights[mask],
	minlength=np.prod(shape)).reshape(shape)

#def mad(a, c=0.6745, axis=0):
	#"""
	#Median Absolute Deviation along given axis of an array:
	#median(abs(a - median(a))) / c
	#"""
	#a = N.asarray(a, N.float64)
	#d = median(a, axis=axis)
	#d = unsqueeze(d, axis, a.shape)
	#return median(N.fabs(a - d) / c, axis=axis)
#import numpy.ma as ma
##from scipy.stats import norm, median
##from scipy.stats.stats import nanmedian,_nanmedian
        
def weighted_avg_and_std(values, weights):
  """
  Return the weighted average and standard deviation.

  values, weights -- Numpy ndarrays with the same shape.
  """
  average = np.average(values, weights=weights)
  variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
  return (average, np.sqrt(variance))

def weighted_percentiles(data, wt, percentiles): 
          """Compute weighted percentiles. 
          If the weights are equal, this is the same as normal percentiles. 
          Elements of the C{data} and C{wt} arrays correspond to 
          each other and must have equal length (unless C{wt} is C{None}). 
   
          @param data: The data. 
          @type data: A L{numpy.ndarray} array or a C{list} of numbers. 
          @param wt: How important is a given piece of data. 
          @type wt: C{None} or a L{numpy.ndarray} array or a C{list} of numbers. 
                  All the weights must be non-negative and the sum must be 
                  greater than zero. 
          @param percentiles: what percentiles to use.  (Not really percentiles, 
                  as the range is 0-1 rather than 0-100.) 
          @type percentiles: a C{list} of numbers between 0 and 1. 
          @rtype: [ C{float}, ... ] 
          @return: the weighted percentiles of the data. 
          """ 
          assert np.greater_equal(percentiles, 0.0).all(), "Percentiles less than zero" 
          assert np.less_equal(percentiles, 1.0).all(), "Percentiles greater than one" 
          data = np.asarray(data) 
          assert len(data.shape) == 1 
          if wt is None: 
                  wt = np.ones(data.shape, np.float) 
          else: 
                  wt = np.asarray(wt, np.float) 
                  assert wt.shape == data.shape 
                  assert np.greater_equal(wt, 0.0).all(), "Not all weights are non-negative." 
          assert len(wt.shape) == 1 
          n = data.shape[0] 
          assert n > 0 
          i = np.argsort(data) 
          sd = np.take(data, i, axis=0) 
          sw = np.take(wt, i, axis=0) 
          aw = np.add.accumulate(sw) 
          if not aw[-1] > 0: 
                  raise ValueError("Nonpositive weight sum")
          w = (aw-0.5*sw)/aw[-1] 
          spots = np.searchsorted(w, percentiles) 
          o = [] 
          for (s, p) in zip(spots, percentiles): 
                  if s == 0: 
                          o.append(sd[0]) 
                  elif s == n: 
                          o.append(sd[n-1]) 
                  else: 
                          f1 = (w[s] - p)/(w[s] - w[s-1]) 
                          f2 = (p - w[s-1])/(w[s] - w[s-1]) 
                          assert f1>=0 and f2>=0 and f1<=1 and f2<=1 
                          assert abs(f1+f2-1.0) < 1e-6 
                          o.append(sd[s-1]*f1 + sd[s]*f2) 
          return o

def weighted_median(data,wt):
  if len(data)==0:
    return np.nan
  return weighted_percentiles(data,wt,[.5])[0]

def weighted_mad(data,wt,c=0.6745):
    a = data
    d = weighted_median(a,wt)
    m = weighted_median(np.abs(a - d) / c,wt)
    return m
  

def mad(a, c=0.6745, axis=None):
    """
    Median Absolute Deviation along given axis of an array:
    median(abs(a - median(a))) / c
    c = 0.6745 is the constant to convert from MAD to std; it is used by
    default
    """
    a = ma.masked_where(a!=a, a)
    if a.ndim == 1:
        d = ma.median(a)
        m = ma.median(ma.fabs(a - d) / c)
    else:
        d = ma.median(a, axis=axis)
        # I don't want the array to change so I have to copy it?
        if axis > 0:
            aswp = ma.swapaxes(a,0,axis)
        else:
            aswp = a
        m = ma.median(ma.fabs(aswp - d) / c, axis=0)
    return m

def nanmedian(arr, **kwargs):
    """
    Returns median ignoring NAN
    """
    return ma.median( ma.masked_where(arr!=arr, arr), **kwargs )

def com(mat,axis=0):
  #TODO
  pass


def COM(Img, xv=[], yv=[], firstOnly=True):
  sz = np.shape(Img)
  if not xv:
    xv = np.arange(sz[1])-sz[1]/2.
  if not yv:
    yv = np.arange(sz[0])-sz[0]/2.

  projX = np.sum(Img,0) 
  projY = np.sum(Img,1)

  # 0th moment ("the mass").
  M = np.sum(Img)

  # 1st moment 
  M1x = np.sum( projX * xv )
  M1y = np.sum( projY * yv )
  # The centre of mass
  COMx = M1x/M;
  COMy = M1y/M;
  # 2nd moment 
  # M2x = sum( projX .* xv.^2 );
  # M2y = sum( projY .* yv.^2 );
  #%% The inertia
  # INERTx = M2x/M;     % Formally, Ix == M2x?
  # INERTy = M2y/M;
  #
  # A variance measure.
  VARx = np.sum( projX * (xv - COMx)**2 ) / M 
  VARy = np.sum( projY * (yv - COMy)**2 ) / M
  SDEVx = np.sqrt(VARx)
  SDEVy = np.sqrt(VARy)

  if firstOnly:
    cm = [M,COMx,COMy,SDEVx,SDEVy]
  else:


    # 3rd moment 
    # A skewness measure.
    # positive --> tail to positive end. zero --> symmetric.
    SKEWx = np.sum( projX * ((xv - COMx)/SDEVx)**3 ) / M 
    SKEWy = np.sum( projY * ((yv - COMy)/SDEVy)**3 ) / M

    # 4th moment 
    # A kurtosis measure.
    # positive --> Matterhorn, negative --> Table Mountain.
    KURTx = np.sum( projX * ((xv - COMx)/SDEVx)**4 ) / M - 3 
    KURTy = np.sum( projY * ((yv - COMy)/SDEVy)**4 ) / M - 3


    cm = [M,COMx,COMy,SDEVx,SDEVy,SKEWx,SKEWy,KURTx,KURTy]
  return cm

def COM1d(I,xv=[],firstOnly=True):
  sz = len(I)
  if len(xv)==0:
    xv = np.arange(sz[1])-sz[1]/2.

  # 0th moment ("the mass").
  M = np.sum(I)

  # 1st moment 
  M1 = np.sum( I * xv )
  # The centre of mass
  COM = M1/M;
  # 2nd moment 
  # M2x = sum( projX .* xv.^2 );
  # M2y = sum( projY .* yv.^2 );
  #%% The inertia
  # INERTx = M2x/M;     % Formally, Ix == M2x?
  # INERTy = M2y/M;
  #
  # A variance measure.
  VAR = np.sum( I * (xv - COM)**2 ) / M 
  SDEV = np.sqrt(VAR)

  if firstOnly:
    cm = [M,COM,SDEV]
  else:


    # 3rd moment 
    # A skewness measure.
    # positive --> tail to positive end. zero --> symmetric.
    SKEW = np.sum( I * ((xv - COM)/SDEV)**3 ) / M 

    # 4th moment 
    # A kurtosis measure.
    # positive --> Matterhorn, negative --> Table Mountain.
    KURT = np.sum( I * ((xv - COM)/SDEV)**4 ) / M - 3 
    cm = [M,COM,SDEV,SKEW,KURT]
  return cm


#def gauss(par=dict(A=[],pos=[],sig=[]),dat=[]):
	#res = par['A']*np.exp(-(dat-par['pos'])**2/par['sig']**2/2)
	#return res

def erfstep(par=dict(h=1,pos=0,sig=1,offs=0),dat=np.linspace(-5,5,100)):
	#res = par['h']*np.exp(-(dat-par['pos'])**2/par['sig']**2/2)
  res = par['h']*(erf((dat-par['pos'])/np.sqrt(2) /par['sig']) + 1)/2 + par['offs']
  return res

def gauss_amp(X,xdat):
  ydat = X[0]*np.exp(-(xdat-X[1])**2/2/X[2]**2)
  return ydat

#def gauss_norm(X,xdat):
#	ydat = 1./sqrt(2.*pi*X[2]**2)*X[0]*exp(-(xdat-X[1])**2/2/X[2]**2)
#	return ydat

def gauss_norm(x,area=1,pos=0,sig=1):
  ydat = 1./np.sqrt(2.*np.pi*sig**2)*area*np.exp(-(x-pos)**2/2/sig**2)
  return ydat
#def gauss_norm(par=dict(area=1,pos=0,sig=1),xdat=np.linspace(-5,5,100)):
  #ydat = 1./np.sqrt(2.*np.pi*par['sig']**2)*par['area']*np.exp(-(xdat-par['pos'])**2/2/par['sig']**2)
  #return ydat


def errortube():
  pass

def histEdges(x,fac=20.):
  lims0 = matplotlib.mlab.prctile(x,p=(20,80))
  ind = (x>lims0[0])&(x<lims0[1])
  interval = np.diff(lims0)/np.round(sum(ind)/fac)
  edges = np.arange(np.min(x),np.max(x),interval)
  return edges

def histogramSmart(x,fac=20.,include=-1,remove=0,maxints=1000000):
  lims0 = np.percentile(x,[20,80])
  ind = (x>lims0[0])&(x<lims0[1])
  interval = np.diff(lims0)/np.round(sum(ind)/fac)
  include = iterfy(include)
  remove = iterfy(remove)
  if sum(include)>0:
    med = np.median(x)
    if len(include)==1:
      include = np.abs(include[0])
      hmn,hmx = np.percentile(x,[50-include,50+include])
    elif len(include)==2:
      hmn,hmx = np.percentile(x,[50-include[0],50+include[1]])
  elif sum(remove)>0:
    med = np.median(x)
    if len(remove)==1:
      remove = np.abs(remove[0])
      hmn,hmx = np.percentile(x,[remove,100-remove])
    elif len(remove)==2:
      hmn,hmx = np.percentile(x,[remove[0],100-remove[1]])
  else:
    hmn,hmx = (np.min(x),np.max(x))

  xd = np.diff(x[ind])
  xd = xd[xd>0]
  xdmn = np.min(xd)
  if xdmn>interval:
    interval = xdmn
    
  if (hmx-hmn)/interval>maxints:
    logbook("Warning: the assigned binwidth %g leads to more bins than assigned in maxint (%g)." %(interval,maxints),level=2)
    interval = (hmx-hmn)/maxints
    logbook("binwidth is set to %g." %(interval))
  edges = np.arange(hmn,hmx,interval)
  h,dum = np.histogram(x,bins=edges)
  return h,edges

def histogram2dSmart(x,y,fac=400,include=-1,remove=0,maxints=500):
  limsx0 = matplotlib.mlab.prctile(x,p=(20,80))
  limsy0 = matplotlib.mlab.prctile(y,p=(20,80))
  indx = (x>limsx0[0])&(x<limsx0[1])
  indy = (y>limsy0[0])&(y<limsy0[1])
  intervalx = np.diff(limsx0)/np.round(sum(indx)/fac)
  intervaly = np.diff(limsy0)/np.round(sum(indy)/fac)
  include = iterfy(include)
  remove = iterfy(remove)
  if sum(include)>0:
    medx = np.median(x)
    medy = np.median(y)
    if len(include)==1:
      includesingle = np.abs(include[0])
      xhmn,xhmx = np.percentile(x,[50-includesingle,50+includesingle])
      yhmn,yhmx = np.percentile(y,[50-includesingle,50+includesingle])
    elif len(include)==2:
      if len(include[0])==1:
        xhmn,xhmx = np.percentile(x,[50-include[0][0],50+include[0][0]])
      elif len(include[0])==2:
        xhmn,xhmx = np.percentile(x,[50-include[0][0],50+include[0][1]])
      if len(include[1])==1:
        yhmn,yhmx = np.percentile(y,[50-include[1][0],50+include[1][0]])
      elif len(include[0])==2:
        yhmn,yhmx = np.percentile(x,[50-include[1][0],50+include[1][1]])
  else:
    xhmn,xhmx = np.min(x),np.max(x)
    yhmn,yhmx = np.min(y),np.max(y)


  if sum(remove)>0:
    medx = np.median(x)
    medy = np.median(y)
    if len(remove)==1:
      removesingle = np.abs(remove[0])
      xhmn,xhmx = np.percentile(x,[removesingle,100-removesingle])
      yhmn,yhmx = np.percentile(y,[removesingle,100-removesingle])
      logbook("here",removesingle)
    elif len(remove)==2:
      remove[0] = iterfy(remove[0])
      remove[1] = iterfy(remove[1])
      if len(remove[0])==1:
        xhmn,xhmx = np.percentile(x,[remove[0][0],100-remove[0][0]])
      elif len(remove[0])==2:
        xhmn,xhmx = np.percentile(x,[remove[0][0],100-remove[0][1]])
      if len(remove[1])==1:
        yhmn,yhmx = np.percentile(y,[remove[1][0],100-remove[1][0]])
      elif len(remove[0])==2:
        yhmn,yhmx = np.percentile(x,[remove[1][0],100-remove[1][1]])
  else:
    xhmn,xhmx = np.min(x),np.max(x)
    yhmn,yhmx = np.min(y),np.max(y)
  #elif sum(remove)>0:
    #med = np.median(x)
    #if len(remove)==1:
      #remove = np.abs(remove[0])
      #hmn,hmx = np.percentile(x,[remove,100-remove])
    #elif len(remove)==2:
      #hmn,hmx = np.percentile(x,[remove[0],100-remove[1]])
  #else:
    #hmn,hmx = (np.min(x),np.max(x))include=-1,remove=0,maxints=1000000
  xd = np.diff(x[indx])
  xd = xd[xd>0]
  xdmn = np.min(xd)
  if xdmn>intervalx:
    intervalx = xdmn
  yd = np.diff(y[indy])
  yd = yd[yd>0]
  ydmn = np.min(yd)
  if ydmn>intervaly:
    intervaly = ydmn
  if (xhmx-xhmn)/intervalx>maxints:
    logbook("Warning: the assigned x binwidth %g leads to more bins than assigned in maxint (%g)." %(intervalx,maxints),level=2)
    intervalx = (xhmx-xhmn)/maxints
    logbook("binwidth is set to %g." %(intervalx),level=2)
  if (yhmx-yhmn)/intervaly>maxints:
    logbook("Warning: the assigned y binwidth %g leads to more bins than assigned in maxint (%g)." %(intervaly,maxints),level=2)
    intervaly = (yhmx-yhmn)/maxints
    logbook("binwidth is set to %g." %(intervaly),level=2)
  edgesx = np.arange(xhmn,xhmx,intervalx)
  edgesy = np.arange(yhmn,yhmx,intervaly)
  h,dumx,dumy = np.histogram2d(x,y,[edgesx,edgesy])
  h = h.transpose()
  return h,edgesx,edgesy

def histVec(v,oversample=1):
  v = np.atleast_1d(v)
  v = np.unique(v)
  vd = np.diff(v)
  vd = np.hstack([vd[0],vd])
  #vv = np.hstack([v-vd/2,v[-1]+vd[-1]/2])
  vv = np.hstack([v-vd/2.,v[-1]+vd[-1]/2.])
  if oversample>1:
    vvo = []
    for i in range(len(vv)-1):
      vvo.append(np.linspace(vv[i],vv[i+1],oversample+1)[:-1])
    vvo.append(vv[-1])
    vv = np.array(np.hstack(vvo))
  return vv

def histVecLinlog(v,smBinsize,lgBinsize):
  V = np.array(v)
  #V.sort()
  dV = np.diff(V)
  smStep = np.min(dV)
  linsteps = dV<(1+.05)*smStep
  indMid = linsteps.nonzero()[-1]+1
  np.arange(V[0]-smBinsize/2.,V[ndMid]+smBinsize,smBinsize)
  

def histVecLinlogFromTT(t,tt,ttOffset=None,binSize=None):
  if ttOffset is None:
    pass
  tT = tt + t + ttOffset

  
  T = np.array(t)
  s = T.argsort()
  Ts = T[s]
  dT = np.diff(Ts)
  smStep = np.min(dT)
  logsteps = (dT>(1+.01)*smStep).nonzero()[0] +1
  linsteps = (dT<=(1+.01)*smStep).nonzero()[0]
  linsteps = np.concatenate([linsteps,[linsteps[-1]+1]])
  linedges = np.arange(np.min(Ts[linsteps])-binSize/2.,np.max(Ts[linsteps])+binSize,binSize)
  

  toLin = tT[s[linsteps]].digitize(linedges)
  toLog = tT[s[logsteps]]
  return ixppy.concatenate([toLin,toLog])





def histVecCenter(v):
  v = np.array(v)*1.
  return v[:-1]+np.diff(v)/2.

#def binStat(M,axis=-1,certaintylevel = 0.95):
#
#
 # 
 ## 
 # 
  #Mmean = []; Mstd = [];
#for jj = 1:length(M(1,:))
    #index = find(~isnan(M(:,jj)));
    #Mmean(jj) = mean(M(index, jj));
    #tMstd  = std(M(index,jj));
    #n = length(index);
    #fac = tinv(certaintylevel,n-1)./sqrt(n);
    #Mstd(jj)  = fac*tMstd;
#end

def getEdges(v,desc):
  assert type(desc) is dict, "Expecting dictionary with named parameters"
  if 'lims' in desc.keys():
    lims = np.sort(desc['lims'])
  elif 'limsPerc' in desc.keys():
    limsperc = desc['limsPerc']
    if not np.iterable(limsperc):
      limsperc = [limsperc,100-limsperc]
    else:
      np.sort(limsperc)
    lims = np.percentile(v,limsperc)
  else:
    lims = [np.min(v),np.max(v)]

  if 'Nbins' in desc.keys():
    Nbins = desc['Nbins']
    if 'ispercbins' in desc.keys():
      ispercbins = desc['ispercbins']
      assert type(ispercbins) is bool, 'key ispercbins needs to be of type bool'
    else:
      ispercbins = False
    if ispercbins:
      edges = np.percentile(v,list(np.linspace(limsperc[0],limsperc[1],Nbins+1)))
    else:
      edges = np.linspace(lims[0],lims[1],Nbins+1)

  elif 'interval' in desc.keys():
    interval = desc['interval']
    edges = np.arange(lims[0],lims[1],interval)

  return edges
