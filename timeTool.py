from __future__ import print_function
import pylab as pl
import os
import numpy as np
import tools
from copy import copy
import ixppy
import time
from toolsLog import logbook

standFilter = None

############ TIMING TOOL #############

def TTextractFromRun(run,opalNo=0,referenceCode=162,refthreshold=True,laseroffCode=67,lmonthreshold=True,outputfile=None,filter='standard'):
  d = dataset(run)
  #TODO
  #d.config.cache_filename = outputfile
  #try:
    #d.config.cache.delete()
  #except:
    #print "no cache file found"
  
  refdat  = getattr(d.eventCode,'code'+str(referenceCode))
  ldat = getattr(d.eventCode,'code'+str(laseroffCode))     
  TTextract(d.__dict__['opal'+str(opalNo)],
            refmon=refdat,lmon=ldat,       
            refthreshold=refthreshold,
            lmonthreshold=lmonthreshold,filter=filter)
  d.save()

def standardfilter():
  if globals()["standFilter"] is not None:
    return globals()["standFilter"]
  else:
    path = os.path.abspath(__file__)
    pathname = os.path.dirname(path)
    weights = np.loadtxt(pathname+'/TTfilt_standard.dat')
    filt = dict()
    filt['weights']=weights
    globals()["standFilter"] = filt
  return filt

def TTextract(det,refmon=None,refthreshold=True,lmon=None,lmonthreshold=True,filter='standard',saveoffs=False):
  ixppy.getProfileLimits(det)
  profiles = ixppy.TTextractProfiles(det,refmon=refmon,refthreshold=refthreshold,lmon=lmon,lmonthreshold=lmonthreshold,Nmax=30,calibcycles=0,append_to_det=False,saveoffs=saveoffs)
  if filter=='create':
    filter = ixppy.TTteachFilter(profiles[0])
  elif filter=='standard':
    filter = TTstandardfilter()
  det._add_saved_datafield('TTfiltsettings',filter)

  ixppy.TTextractProfiles(det,refmon=refmon,refthreshold=refthreshold,lmon=lmon,lmonthreshold=lmonthreshold,saveoffs=saveoffs)
  ixppy.TTextractFilterPositions(det)

def TTteachFilter(profiles,):
  from scipy.linalg import toeplitz
  nfh = tools.nfigure('Digital filter design: Select signal')
  pl.clf()
  p = profiles[2,:]
  pl.plot(p)
  siglim = np.round(tools.getSpanCoordinates('horizontal'))
  #print 'Select step position of step'
  #stepoffs = pl.ginput(1)[0][0]-np.mean(siglim)
  #pl.axvline(stepoffs+np.mean(siglim),color='r')
  ys = p[siglim[0]:siglim[1]]
  ysacl = np.correlate(ys,ys,'same')
  sacl = ysacl

  nfh = tools.nfigure('Digital filter design: Select noise region')
  pl.clf()
  pl.plot(profiles.transpose())
  logbook("select lower limit of noise area (NB: has same width as signal range!)")
  noiselim = pl.ginput(1)
  noiselim = round(noiselim[0][0])+np.array([0,np.diff(np.array(siglim))[0]])
  pl.axvspan(noiselim[0],noiselim[1],facecolor='r',alpha=0.5)
  pl.axvline(noiselim[0],color='r')
  pl.axvline(noiselim[1],color='r')
  logbook(noiselim)
  nacl = []
  for p in profiles:
    yn = p[noiselim[0]:noiselim[1]]
    ynacl = np.correlate(yn,yn,'same')
    nacl.append(ynacl)
  nacl = np.mean(np.vstack(nacl),axis=0)

  Ynacl = toeplitz(nacl,r=np.zeros(len(nacl)))
  R  = np.matrix(Ynacl).I
  Rs = np.matrix(sacl)

  weights = R*Rs.transpose()
  weights = np.array(weights).ravel()
  weights = weights-np.median(weights)


  filtsettings = dict(weights=np.array(weights),
                      #stepoffs=stepoffs,
                      noise_limits=noiselim)
  return filtsettings

def applyFilter(data,filtsettings,plotOutput=False,polysettings=None,erfsettings=None,saveplots=False,kind="stepUp"):
  weights = np.array(filtsettings['weights']).ravel()
  #stepoffs = filtsettings['stepoffs'] 
  lf = len(weights)
  halfrange = round(lf/10)
  pos = []
  amp = []
  fwhm = []
  runningno = 0
  for d in data:
    f0 = np.convolve(np.array(weights).ravel(),d,'same')
    f = f0[lf/2:len(f0)-lf/2-1]
    if (kind=="stepUp"):
      mpr = f.argmax()
    else:
      mpr = f.argmin()
    # now do a parabolic fit around the max
    xd = np.arange(max(0,mpr-halfrange),min(mpr+halfrange,len(f)-1))
    yd = f[max(0,mpr-halfrange):min(mpr+halfrange,len(f)-1)]
    p2 = np.polyfit(xd,yd,2)
    tpos = -p2[1]/2./p2[0]
    tamp = np.polyval(p2,tpos)
    try:
      beloh = (f<tamp/2).nonzero()[0]-mpr
      tfwhm = abs(beloh[beloh<0][-1]-beloh[beloh>0][0])
    except:
      logbook("FWHM not applied",level=0)
      tfwhm = np.nan
    pos.append(tpos)
    amp.append(tamp)
    fwhm.append(tfwhm)
    runningno+=1
  pos  = np.asarray(pos) + lf/2.
  amp  = np.asarray(amp)
  fwhm = np.asarray(fwhm)
  returntuple = [pos,amp,fwhm]
  return tuple(returntuple)

def TTextractFilterPositions(Areadet,filtsettings=None,polysettings=None):
  if not filtsettings:
    filtsettings = Areadet.TTfiltsettings
  #if not polysettings:
    #polysettings = dict(rpts=100,cpts=20)
  pos = []
  amp = []
  fwhm = []
  if polysettings:
    poly_pos = []
    poly_cen = []
  ccN = 0
  for ds in Areadet.TTtraces:
    logbook('...extracting from cc %d'%(ccN))
    data = ds[:]
    if polysettings:
      tpos,tamp,tppos,tpcen = TTapplyFilter(data,filtsettings,polysettings=polysettings)
    else:
      tpos,tamp,tfwhm = TTapplyFilter(data,filtsettings)
    pos.append(tpos)
    amp.append(tamp)
    fwhm.append(tfwhm)
    if polysettings:
      poly_pos.append(tppos)
      poly_cen.append(tpcen)
    ccN+=1
  Areadet._add_saved_datafield('TTfiltPos',pos)
  Areadet._add_saved_datafield('TTfiltAmp',amp)
  Areadet._add_saved_datafield('TTfiltFwhm',fwhm)
  if polysettings:
    Areadet._add_saved_datafield('TTfiltPolyPos',poly_pos)
    Areadet._add_saved_datafield('TTfiltPolyCen',poly_cen)



def extractProfiles(Areadet, xrayoff=None, laseroff=None,
    profileLimits=None, transpose=False, saveoffs=False,dataset_name_traces='TTtraces'):
  """
    Areadet is the dataset with the camera images
    'refmon' is the incoming intensity monitor (to check for x-ray off)
    'refthreshold' threshold to find x-ray off
    'Nxoff' number of x-ray off images to average around the image to analyze
    'profileLimits' ROI {'limits': np.array([ 146.,  198.]), 'projection': 'vertical range'} os the kind
    'profileLimitsRef' same as above for the reference trace
    'Nmax' limits to Nmax images per calibcycle
    'steps' analyze only every steps images
    'calibcycles' which ones to do, if None do all
  """
  if not laseroff is None:
    dat = Areadet.data * laseroff.filter(False).ones()
  else:
    dat = Areadet.data
  #dat = dat - 32
  
  proflimits,profile = ixppy.getProfileLimits(dat,lims=profileLimits,transpose=transpose)
  
  datpump  = profile * xrayoff.filter(False).ones()
  datref   = profile * xrayoff.filter(True).ones()
  datrefIP = datref.interpolate(xrayoff.filter(False).time)
  
  Areadet[dataset_name_traces] = (datpump-datrefIP)/datrefIP
  #Areadet['prof'] = profile
  #Areadet['datpump'] = datpump
  #Areadet['datrefIP'] = datrefIP

def extractProfilesCorr(Areadet, xrayoff=None, laseroff=None, dataset_name_traces='TTtraces',dataset_name_traces_raw='TTtraces_raw',evaluate=False):
  """
    Areadet is the dataset with the camera images
    'refmon' is the incoming intensity monitor (to check for x-ray off)
    'refthreshold' threshold to find x-ray off
    'Nxoff' number of x-ray off images to average around the image to analyze
    'profileLimits' ROI {'limits': np.array([ 146.,  198.]), 'projection': 'vertical range'} os the kind
    'profileLimitsRef' same as above for the reference trace
    'Nmax' limits to Nmax images per calibcycle
    'steps' analyze only every steps images
    'calibcycles' which ones to do, if None do all
  """
  if not laseroff is None:
    dat = Areadet[dataset_name_traces_raw] * laseroff.filter(False).ones()
  else:
    dat = Areadet[dataset_name_traces_raw]
  #dat = dat - 32
  
  
  datpump  = dat * xrayoff.filter(False).ones()
  datref   = dat * xrayoff.filter(True).ones()
  datrefIP = datref.interpolate(xrayoff.filter(False).time)
  datrefIP[0,0]
  Areadet[dataset_name_traces] = (datpump-datrefIP)/datrefIP
  if evaluate:
    Areadet[dataset_name_traces].evaluate(force=True)

def extractProfilesOnly(Areadet, profileLimits=None, transpose=False,dataset_name_traces='TTtraces_raw',force=True):
  """ 
    Areadet is the dataset with the camera images
    'refmon' is the incoming intensity monitor (to check for x-ray off)
    'refthreshold' threshold to find x-ray off
    'Nxoff' number of x-ray off images to average around the image to analyze
    'profileLimits' ROI {'limits': np.array([ 146.,  198.]), 'projection': 'vertical range'} os the kind
    'profileLimitsRef' same as above for the reference trace
    'Nmax' limits to Nmax images per calibcycle
    'steps' analyze only every steps images
    'calibcycles' which ones to do, if None do all
  """
  #if not laseroff is None:
    #dat = Areadet.data * laseroff.filter(False).ones()
  #else:
  dat = Areadet.data
  #dat = dat - 32
  
  proflimits,profile = ixppy.getProfileLimits(dat,lims=profileLimits,transpose=transpose)
  
  Areadet[dataset_name_traces] = profile
  Areadet[dataset_name_traces].evaluate(force=force)

def extractFromRunList(runlist,exp,datasetname='opal2',profileLimits=None,xrayoffCode=None,laseroffCode=None,filter=None,save=False):
  for run in runlist:
    d = ixppy.dataset((exp,run))
    logbook("TT extracting from run %d" %run)
    extractFromRun(d,datasetname=datasetname,profileLimits=profileLimits,xrayoffCode=xrayoffCode,laseroffCode=laseroffCode,filter=filter,save=save)
    logbook("done!")

    

def extractFromRun(d,datasetname='opal2',profileLimits=None,xrayoffCode=None,laseroffCode=None,filter=None,save=False):
  dataset = d[datasetname]
  extractProfilesOnly(dataset,profileLimits)
  xrayoff = d.eventCode['code_%s'%xrayoffCode]
  if laseroffCode is not None:
    laseroff = d.eventCode['code_%s'%laseroffCode]
  else:
    laseroff = None

  extractProfilesCorr(dataset,xrayoff=xrayoff,laseroff=laseroff,evaluate=True)
  applyFilterToAll(dataset)
  if save:
    d.save()

def applyFilterToAll(det,filter=None,kind='stepUp'):
  if filter is None:
    filter = standardfilter()
  traces = det.TTtraces
  o = ixppy.applyFunction(applyFilter,[traces,filter],dict(kind=kind),outputtypes=['memdata']*3,forceCalculation=True,transposeStack=False)
  det['TTpos'] = o[0]
  det['TTamp'] = o[1]
  det['TTfwhm'] = o[2]


  #Areadet[dataset_name_traces] = dat




  #d.save()

  # Plotting result
  #allcctraces_stacked = np.vstack(allcctraces) 
  #pl.imshow(allcctraces_stacked,interpolation='nearest')
  #pl.axis('normal');pl.axis('tight')
  #pl.draw()

def TTcalc_weightedRatio(det,mon,TTdet,tvec=None):
  #if not tvec:
    #tvec = 
  timevec = []
  for tvecS in tvec:
    timevec.append(tvecS + 1e-12*np.polyval(TTdet.TTpxCalib,TTdet.TTfiltPos))
  wR = calc_weightedRatio(timevec,det,mon)

#############

