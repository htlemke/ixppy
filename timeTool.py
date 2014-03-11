import pylab as pl
#import h5py
import os
#from os.path import expanduser
#from socket import gethostname
#import dateutil
#import sys
import numpy as np
#import ixppy_specialdet
import tools
#from toolsVarious import addToObj
#import toolsHdf5 as tH5
#from toolsHdf5 import datasetRead as h5r
#from functools import partial,wraps
from copy import copy
#import re
#import examples
import ixppy
#import datetime
#import operator
import time
#import lclsH5
#import copy as pycopy


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
  print "select lower limit of noise area (NB: has same width as signal range!)"
  noiselim = pl.ginput(1)
  noiselim = round(noiselim[0][0])+np.array([0,np.diff(np.array(siglim))[0]])
  pl.axvspan(noiselim[0],noiselim[1],facecolor='r',alpha=0.5)
  pl.axvline(noiselim[0],color='r')
  pl.axvline(noiselim[1],color='r')
  print noiselim
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
  if polysettings:
    poly_pos = []
    poly_cen = []
  for d in data:
    #print runningno
    #runningno+=1
    try:
      if max(abs(d))>5:
        raise Exception("Strange array!")
      f0 = np.convolve(np.array(weights).ravel(),d,'same')
      f = f0[lf/2:len(f0)-lf/2-1]
      if (kind=="stepUp"):
        mpr = f.argmax()
      else:
        mpr = f.argmin()
    #if True:
      xd = np.arange(max(0,mpr-halfrange),min(mpr+halfrange,len(f)-1))
      yd = f[max(0,mpr-halfrange):min(mpr+halfrange,len(f)-1)]
      p2 = np.polyfit(xd,yd,2)
      
      tpos = -p2[1]/2./p2[0]
      tamp = np.polyval(p2,tpos)
      try:
        beloh = (f<tamp/2).nonzero()[0]-mpr
        tfwhm = abs(beloh[beloh<0][-1]-beloh[beloh>0][0])
      except:
        print "FWHM not applied"
        tfwhm = np.nan

      #tools.nfigure('test')
      #pl.subplot(211)
      #pl.plot(d)
      #pl.hold(True)
      #pl.axvline(tpos)
      ##pl.axvline(tpos+lf/2+stepoffs,color='r')
      #pl.axvline(tpos+lf/2,color='r')
      #pl.axvline(tpos+lf/2)
      #pl.axvline(tpos+lf)
      #pl.hold(False)
      #pl.subplot(212)
      #pl.plot(f0)
      #pl.hold(True)
      #pl.axvline(tpos)
      ##pl.axvline(tpos+lf/2+stepoffs,color='r')
      #pl.axvline(tpos+lf/2,color='r')
      #pl.axvline(tpos+lf/2)
      #pl.axvline(tpos+lf)
      #pl.hold(False)
      #pl.draw()
      #pl.waitforbuttonpress()

      if polysettings:
        tcen = tpos+lf/2
        rpts = polysettings['rpts']
        cpts = polysettings['cpts']
        # make region for polyfit
        pxd = np.arange(int(tcen-rpts),int(tcen+rpts))
        pyd = d[int(tcen-rpts):int(tcen+rpts)]
        #import pdb;pdb.set_trace() 
        # do the fit to find step on first order
        p10 = np.polyfit(pxd,pyd,10)
        dp10= np.polyder(p10)
        cpxd = pxd[cpts:-cpts]
        cpyd = pyd[cpts:-cpts]
        pkpos0_ind = np.polyval(dp10,cpxd).argmax()

        # go for extrema
        ddp10 = np.polyder(dp10)
        
        zs = np.roots(ddp10) - cpxd[pkpos0_ind]
        zs = zs[np.imag(zs)==0]
        zs = np.real(zs)
        pkpos_ind = np.argmin(np.abs(zs))
        pkpos = zs[pkpos_ind]  + cpxd[pkpos0_ind]
        #limits for erf fit
        #print np.roots(dp10) -cpxd[pkpos0_ind]
        zs = np.roots(dp10) - cpxd[pkpos0_ind]
        zs = zs[np.imag(zs)==0]
        zs = np.real(zs)
        rp = np.min(zs[zs>0]) + cpxd[pkpos0_ind]
        lp = np.max(zs[zs<0]) + cpxd[pkpos0_ind]
        tpoly_pos = pkpos
        tpoly_cen = (rp+lp)/2.
        #print lp,rp
        if False:
          steprad = (rp-lp)/2
          efx = np.float64( np.arange(int(lp-steprad),int(rp+steprad)))
          efy = d[int(lp-steprad):int(rp+steprad)]
          #de=bug
          startpar = [np.mean(efy[-3:])-np.mean(efy[:3]),
                                   np.mean(efx),
                                   (rp-lp)/10.,
                                   np.mean(efy[:3])]
          mfh = tools.minuitfit(tools.erfstep,startpar,efx,efy)
          mfh.migrad()
          #erffitres = mfh.values
        # find maximum on finer procedure (should be not necessary in a while)

        #ccpxd = cpxd[pkpos0_ind-5:pkpos0_ind+5]
        #ccpyd = cpyd[pkpos0_ind-5:pkpos0_ind+5]
        #pp2 = np.polyfit(ccpxd,ccpyd,2)
        #tpoly_pos = -pp2[1]/2./pp2[0]
        #poly_pos.append(tpoly_pos)
         
        ### PLOTTING
        print plotOutput
        if plotOutput:
          tools.nfigure('test')
          mah = pl.subplot(211)
          pl.plot(d)
          pl.hold(True)
          pl.plot(pxd,np.polyval(p10,pxd),'r--')
          pl.plot(pxd,
                  np.polyval(dp10,pxd)/np.max(np.abs(np.polyval(dp10,pxd)))*.04,
                  'm')
          pl.axvline(tpos,ls='--')
          #pl.axvline(tpos+lf/2+stepoffs,color='r')
          pl.axvline(tpos+lf/2,ls='--')
          pl.axhline(0,color='k')

          pl.axvline(tpos+lf,ls='--')
          pl.axvline(pkpos,color='c')
          pl.axvline(lp,color='y')
          pl.axvline(rp,color='y')
          #pl.plot(efx,tools.erfstep(par=dict(h=startpar[0],
                                      #pos=startpar[1],
                                      #sig=startpar[2],
                                      #offs=startpar[3]),dat=efx),'r')
          #pl.plot(efx,efy,'m')
          pl.hold(False)
          
          pl.subplot(212,sharex = mah)
          pl.plot(f0)
          pl.hold(True)
          pl.axvline(tpos)
          #pl.axvline(tpos+lf/2+stepoffs,color='r')
          pl.axvline(tpos+lf/2)
          pl.axvline(pkpos,color='c')
          pl.axvline(lp,color='y')
          pl.axvline(rp,color='y')
          pl.hold(False)
          pl.draw()



        ### PLOTTING
      
      if plotOutput:
        tools.nfigure('filter')
        mah = pl.subplot(211)
        pl.plot(d)
        pl.hold(True)
        pl.axvline(tpos)
        #pl.axvline(tpos+lf/2+stepoffs,color='r')
        pl.axvline(tpos+lf/2,color='r')
        pl.axvline(tpos+lf/2)
        pl.axvline(tpos+lf)
        pl.hold(False)
	pl.ylabel('Rel. transmission change')
        pl.subplot(212,sharex = mah)

        yf = np.polyval(p2,xd)
        pl.plot(f,'k')
        pl.hold(True)
        pl.plot(xd,yf,'r')
        pl.axvline(mpr,color='k')
        pl.axvline(tpos,color='r')
        pl.axhline(tamp,color='r')
        pl.axhline(tamp/2,color='k')
        pl.axhline(0,color='k')
        pl.axvline(mpr-lf/2,color='g')
        pl.axvline(mpr+lf/2,color='g')
        pl.axvline(mpr-lf/4,color='c')
        pl.axvline(mpr+lf/4,color='c')
        if not np.isnan(tfwhm):
          pl.axvline(beloh[beloh<0][-1]+mpr,color='b')
          pl.axvline(beloh[beloh>0][0]+mpr,color='b')
        if polysettings:
          pl.axvline(pkpos-lf/2,color='y')
          pl.axvline(rp-lf/2,color='y')
          pl.axvline(lp-lf/2,color='y')

        pl.hold(False)
	pl.ylabel('Filtered')
	pl.xlabel('Spectral bin / px')
        pl.draw()
	if saveplots:
	  pl.gcf().savefig('%s_%04d.png'%(saveplots,runningno))
	else:
          pl.waitforbuttonpress()

    except Exception, e:
      print e
    #else:
      tpos = 0
      tamp = 0
      tfwhm = 0
      if polysettings:
        tpoly_pos = 0
        tpoly_cen = 0

    pos.append(tpos)
    amp.append(tamp)
    fwhm.append(tfwhm)
    if polysettings:
      poly_pos.append(tpoly_pos)
      poly_cen.append(tpoly_cen)
    runningno+=1

  pos = np.array(pos)
  amp = np.array(amp)
  fwhm = np.array(fwhm)
  returntuple = [pos,amp,fwhm]
  if polysettings:
    poly_pos = np.array(poly_pos)
    poly_cen = np.array(poly_cen)
    returntuple.append(poly_pos)
    returntuple.append(poly_cen)

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
    print '...extracting from cc %d'%(ccN)
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
  if not laseroff == None:
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

def extractProfilesCorr(Areadet, xrayoff=None, laseroff=None, dataset_name_traces='TTtraces',dataset_name_traces_raw='TTtraces_raw'):
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
  if not laseroff == None:
    dat = Areadet[dataset_name_traces_raw] * laseroff.filter(False).ones()
  else:
    dat = Areadet[dataset_name_traces_raw]
  #dat = dat - 32
  
  
  datpump  = dat * xrayoff.filter(False).ones()
  datref   = dat * xrayoff.filter(True).ones()
  datrefIP = datref.interpolate(xrayoff.filter(False).time)
  
  Areadet[dataset_name_traces] = (datpump-datrefIP)/datrefIP

def extractProfilesOnly(Areadet, profileLimits=None, transpose=False,dataset_name_traces='TTtraces_raw'):
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
  #if not laseroff == None:
    #dat = Areadet.data * laseroff.filter(False).ones()
  #else:
  dat = Areadet.data
  #dat = dat - 32
  
  proflimits,profile = ixppy.getProfileLimits(dat,lims=profileLimits,transpose=transpose)
  
  Areadet[dataset_name_traces] = profile
  Areadet[dataset_name_traces].evaluate()


def applyFilterToAll(det,filter=None):
  if filter==None:
    filter = standardfilter()
  traces = det.TTtraces
  o = ixppy.applyFunction(applyFilter,[traces,filter],dict(),outputtypes=['memdata']*3,forceCalculation=True)
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

