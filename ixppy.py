import pylab as pl
import h5py
import os
from os.path import expanduser
from socket import gethostname
import dateutil
import sys
import numpy as np
#import ixppy_specialdet
import tools
from toolsVarious import addToObj
import toolsHdf5 as tH5
from toolsHdf5 import datasetRead as h5r
from functools import partial,wraps
from copy import copy
import re
import examples
import ixppy
import datetime
import operator
import time
import lclsH5
import copy as pycopy
try:
  import psana
except:
  psana = None
  print "psana not available on present machine, hdf5 files required."

import progressbar as pb

############## DATASET ############################

class dataset(object):
  def __init__(self,
    inputFilesOrExpRunTuple='',
    detectors = [],
    ixpFile=None,
    beamline = None,               
    rdPointDetectorsImmediately=True,
    rdTimestamps=True,
    sortForTimestamps=True,
    outputDirectory = '',
    readCachedData = True,
    cacheFile = None,
    ignore_daqHdf5file = False
              ):
    self._name = 'dataset'
    # dropObject is a 'container'
    self.config = tools.dropObject()
    # read configuration file
    self._rdConfiguration(beamline)
    # boolean if cache data is taken into account 
    # HDF5/XTC/IXP file(s) to read, it is a tuple of lists for files in the three formats, empty if files not available.
    allfilenames   = self._getFilename(inputFilesOrExpRunTuple)
    self.config.fileNamesH5 = allfilenames[0]
    self.config.fileNamesXtc = allfilenames[1]
    self.config.fileNamesIxp = allfilenames[2]
    # check available formats, decide on strategy to use, user input possible, sets self.config.filestrategy variable.
    self.config.readCachedData = readCachedData
    self._getFileStrategy()
    if len(self.config.fileNamesIxp)==0:
      if ixpFile==None:
	tfname =  self.config.fileNamesH5[0]
	ixpname = os.path.join(self.config.cachePath,os.path.basename(tfname))
	ixpname = os.path.splitext(ixpname)[0] + '.ixp.h5'
	self.config.ixp = Ixp(ixpname)
      elif os.path.isdir(ixpFile):
	tfname =  self.config.fileNamesH5[0]
	ixpname = os.path.join(ixpFile,os.path.basename(tfname))
	ixpname = os.path.splitext(ixpname)[0] + '.ixp.h5'
	self.config.ixp = Ixp(ixpname)

      else:
	self.config.ixp = Ixp(ixpFile)
	if os.path.exists(ixpFile):
	    self.config.filestrategy.append('ixp')
    else:
      if ixpFile==None:
        self.config.ixp = Ixp(self.config.fileNamesIxp[0])
      else:
        self.config.ixp = Ixp(ixpFile)
	if os.path.exists(ixpFile):
	  self.config.filestrategy.append('ixp')

    self._ixpHandle = self.config.ixp
    self._ixpsaved = []
    if 'ixp' in self.config.filestrategy:
      self.config.ixp.get_cacheFileHandle(reopen=True)
      dat = self.config.ixp.load()
      if 'dataset' in [ixps[0] for ixps in dat._ixpsaved]:
	names = [tn[0] for tn in dat.dataset._ixpsaved]
	for name in names:
	  self.__setitem__(name,dat.dataset[name],setParent=True)

    
    if 'h5' in self.config.filestrategy:
      self.config.lclsH5obj= lclsH5.lclsH5(self.config.fileNamesH5,self.config.cnfFile)
      
      self.config.lclsH5obj.checkFiles()
      self.config.lclsH5obj.findDetectors(detectors)
      self.config.lclsH5obj.initDetectors()
      if hasattr(self.config.lclsH5obj,'scanVars'):
	for name in self.config.lclsH5obj.scanVars.keys():
	  tVar = self.config.lclsH5obj.scanVars[name]
	  if not hasattr(tVar,'names'): continue
	  for vname,vdat in zip(tVar.names,np.asarray(tVar.data).T):
	    vname = getValidFieldname(vname,lowerit=True) 
            addToObj(self,name+'.'+vname,vdat,ixpsaved=True)
      if hasattr(self,'scan'):
	scan=self.scan
      else:
	scan=None

      for detName in self.config.lclsH5obj.detectorsNames:
        det = self.config.lclsH5obj.detectors[detName]
        if det._isPointDet:
	  if hasattr(det,'fields'):
	    for fieldName in det.fields.keys():
              addToObj(self,detName+'.'+fieldName,memdata(name=fieldName,input = det.fields[fieldName],scan=scan),ixpsaved=True)
	else:
	  #self[detName] = tools.dropObject(parent=self)
	  #self[detName].data = data(name=detName,time=det.time,input = det.readData,parent=self[detName],scan=scan)
          addToObj(self,detName+'.data',data(name=detName,time=det.time,input = det.readData,scan=scan),ixpsaved=True,setParent=True)

	  

  def _add(self,name,data,ixpsaved='auto'):
    self.__dict__[name]=data
    self._ixpsaved.append((name,ixpsaved))
  #def __repr__(self):
    #return "dropObject with fields: "+str(self.__dict__.keys())
  def __getitem__(self,x):
    return self.__dict__[x]
  def __setitem__(self,name,var,setParent=True):
    self._add(name,var)
    if setParent:
      try:
        self[name]._parent = self
      except:
	pass

  def _getItemsNameList(self,excludeFun=True,excludeConf=True,excludePrivate=True):
    lst = self.__dict__.keys()
    if excludePrivate:
      for l in lst:
	if l[0]=='_':
	  lst.remove(l)
    #if excludeFun:
      #for l in lst:
	#if hasattr('__call__',self.):
	  #lst.remove(l)
    if excludeConf:
      if 'config' in lst:
	lst.remove('config')
    return lst

  def items(self):
    names = self._getItemsNameList()
    op = []
    for name in names:
      op.append((name,self[name]))
    return op
  
  def itemNames(self):
    return self._getItemsNameList()
	  
  def save(self,name=None,force=False):
    #self.config.ixp.save(self,self._ixpHandle,name='dataset')
    self._ixpHandle.save(self,self._ixpHandle.fileHandle,name='dataset',force=force)
    #self._checkCalibConsistency()
  # END OF DATASET.__init__

  def _checkCalibConsistency(self):
    import itertools
    Nc = []
    # upon initialization each detector check for Ncalib
    tocheck = list(itertools.chain( self.detectors.values(),self._scanVars))
    Nc = []
    for d in tocheck:
      Nc.append(d._numOfScanSteps)
    Nc = np.asarray(Nc)
    NcMin = Nc.min(axis=0)
    NcMax = Nc.max(axis=0)
    for d in tocheck:
      # print out most limiting detector
      if (list(NcMin) == list(d._numOfScanSteps)) and (list(NcMin)!=list(NcMax)):
        print "WARNING: Detector/Scan ",d,"is limiting the number of Calybcycle to",str(NcMin),"instead of ",str(NcMax)
      d._numOfScanSteps = list(NcMin)
    self.numOfScanSteps = list(NcMin)
    if len(NcMin) ==1:
      self.numOfScanSteps = self.numOfScanSteps[0]

  def _findData(self,subSelection=[]):
    """ Finds detectors in hdf5 file matching with mnemonic given in config file;
    the matching mnemonic names are as dictionaries (self.pointDet and self.areaDet)
    The 
    """
    if (subSelection==[]) or (subSelection is None):
      subSelection = self.config.cnfFile["pointDet"].keys() + self.config.cnfFile["areaDet"].keys()
    h = self.config.fileHandlesH5[0]
    pointDet = self.config.cnfFile["pointDet"]
    # try to use only CalibCycle0
    try:
      base = "Configure:0000/Run:0000/CalibCycle:0000/"
      h5names = tH5.getDataset(h[base])
      h5names = [base+x for x in h5names]
      # find all confs
      base = "Configure:0000/"
      confs = h[base].keys()
      h5confs = []
      for c in confs:
        if (c.find("Run")==0):
          continue
        else:
          temp = tH5.getDataset(h[base][c])
          for t in temp:
            h5confs.append(base+c+"/"+t)
    except KeyError:
      h5names = tH5.getDataset(h)
    ret = {}
    # *** start EpicsPV *** #
    # look for epics name
    epicsFound=False
    if ("epics_dset" in self.config.cnfFile):
      epicsMne = self.config.cnfFile["epics_dset"][0]
      epicsReg = self.config.cnfFile["epics_dset"][1]
      epicsH5Names=[x for x in h5names if (x.find(epicsReg)>-1)]
      # common Epics path:
      ntemp = min([len(x.split("/")) for x in epicsH5Names])
      epicsCommon = "/".join(epicsH5Names[0].split("/")[0:ntemp])
      # epics var
      self._epicsPaths = {}
      for d in h[epicsCommon]:
        dpath = d
        d = d.replace(':','_')
        d = d.replace('-','_')
        d = d.replace(' ','_')
        d = d.replace('.','_')
        mne = "%s.%s" % (epicsMne.split("/")[0],d)
        self._epicsPaths[mne]={}
        self._epicsPaths[mne]["data"] = epicsCommon.replace('CalibCycle:0000','CalibCycle:%04d')+"/"+dpath+"/data"
        self._epicsPaths[mne]["time"] = epicsCommon.replace('CalibCycle:0000','CalibCycle:%04d')+"/"+dpath+"/time"
        self._epicsPaths[mne]["conf"] = []
      self._epicsNames = self._epicsPaths.keys()
    else:
      self._epicsNames = []
    # *** stop EpicsPV *** #
    for (mnemonic,name) in pointDet.iteritems():
      if (mnemonic.find("epics")>-1) and (mnemonic.find("*")>-1):
        continue
      mnemonic = mnemonic.split('_bak')[0]
      # skip if not in the group we want to read
      if mnemonic not in subSelection:
        continue
      nameData = name["data"].replace("*","\S+")
      detDataset = [x for x in h5names if (re.search(nameData,x) is not None)]
      nameConf = name["conf"].replace("*","\S+")
      detConf    = [x for x in h5confs if (re.search(nameConf,x) is not None)]
      data = [x for x in detDataset if x[-5:]=="/data"]
      time = [x for x in detDataset if x[-5:]=="/time"]
      if ( (len(data) != 0) and (len(time) != 0) ):
        ret[mnemonic] = {}
        ret[mnemonic]["data"] = data[0].replace('CalibCycle:0000','CalibCycle:%04d')
        ret[mnemonic]["time"] = time[0].replace('CalibCycle:0000','CalibCycle:%04d')
        if len(detConf)>0:
          ret[mnemonic]["conf"] = detConf[0]
    self._pointDetPaths = ret
    self.pointDetNames = ret.keys()
    areaDet = self.config.cnfFile["areaDet"]
    ret = {}
    # 3D detectors need special care because data are written differently 
    # /data, /image, /waveform
    for (mnemonic,name) in areaDet.iteritems():
      mnemonic = mnemonic.split('_bak')[0]
      # skip if not in the group we want to read
      if mnemonic not in subSelection:
        continue
      name = name["data"].replace("*","\S+")
      name_nodata = "/".join(name.split("/")[0:-1])
      detDataset = [x for x in h5names if (re.search(name_nodata,x) is not None)]
      conf = [ ]
      data = [x for x in detDataset if (re.search(name,x) is not None)]
      time = [x for x in detDataset if x[-5:]=="/time"]
      if ( (len(data) != 0) and (len(time) !=0) ):
        ret[mnemonic] = {}
        ret[mnemonic]["data"] = data[0].replace('CalibCycle:0000','CalibCycle:%04d')
        ret[mnemonic]["time"] = time[0].replace('CalibCycle:0000','CalibCycle:%04d')
        ret[mnemonic]["conf"] = conf
    self._areaDetPaths = ret
    self.areaDetNames = ret.keys()
    self._detectorsPaths = tools.dictMerge(self._pointDetPaths,self._areaDetPaths)
    self.detectorsNames = self.pointDetNames + self.areaDetNames
    # *** start scan variables *** #
    temp = []
    if (len(self.config.cnfFile["scan_step"])>0):
      for scan_var in self.config.cnfFile["scan_step"]:
        mne,reg = scan_var
        reg  = reg.replace("*","\S+")
        data = [x for x in h5names if (re.search(reg,x) is not None)]
        path = data[0].replace('CalibCycle:0000','CalibCycle:%04d')
        try:
          obj = scanVar(self.config.fileHandlesH5,mne,path)
          tools.addToObj(self,mne,obj)
          temp.append(obj)
        except:
          pass
    self._scanVars = temp
    # *** stop scan variables *** #
    return

  def _rdConfiguration(self,beamline):
    if not beamline==None:
      self.config.beamline = beamline
      self.config.cnfFile = rdConfiguration(beamline=beamline)
    else:
      self.config.cnfFile = rdConfiguration()
      self.config.beamline = self.config.cnfFile['beamline']
    # setup path for data and cached files
    self.config.hostname = gethostname()
    knownhosts = self.config.cnfFile['dataPath'].keys()
    if self.config.hostname in knownhosts:
       self.config.dataPath = self.config.cnfFile['dataPath'][self.config.hostname]
    else:
       self.config.dataPath = os.path.join(self.config.cnfFile['dataPath']['default'],self.config.beamline)
    knownhosts = self.config.cnfFile['cachePath'].keys()
    if self.config.hostname in knownhosts:
       self.config.cachePath = self.config.cnfFile['cachePath'][self.config.hostname]
    else:
       self.config.cachePath = os.path.join(self.config.cnfFile['cachePath']['default'])

  def _getFilename(self,inputFilesOrExpRunTuple):
    # use a shorter local name
    f = inputFilesOrExpRunTuple
    # interprets input (either filename or
    # check if all string
    StrOrStrList = np.alltrue(tools.iterate(f,isinstance,str))
    if (StrOrStrList):
      filenames = f
    elif type(f) is tuple:
      if type(f[0]) is str:
        """ you are just passing the experiment name, e.g. 'xpp66613' as first arg of f"""
        self.config.experiment = f[0]
      if type(f[0]) is int:
        """ you are just passing the experiment number, e.g. 66613 as first arg of f"""
        self.config.experiment = self.config.beamline+str(f[0])
        
      self.config.run = f[1]
      filenames = self._makeFilenameFromExpRun()
    filenamesH5 =tools.iterfy(filenames)

    # TODO: replace "isavailable" by properties checking for empty list
    for f in filenamesH5:
      if (not os.path.exists(f)):
        print "Asked to read file %s, but is does not exist" % f
        self.config.daqHdf5file_available = False
      else:
        self.config.daqHdf5file_available = True
    
    # TODO get according xtc path and check if present.
    filenamesXtc = []
    
    # check if cached files are present
    if len(filenamesH5)==1 and filenamesH5[0][-7:]=='.ixp.h5':
      filenamesIxp = filenamesH5
      filenamesH5 = []
    else:

      filenamesIxp = []
      for f in filenamesH5:
	cached_filename = os.path.join(self.config.cachePath,os.path.basename(f))
	cached_filename = os.path.splitext(cached_filename)[0] + '.ixp.h5'
	if tools.fileExists(cached_filename):
	  filenamesIxp.append(cached_filename)
      
    return (filenamesH5,filenamesXtc,filenamesIxp)


  def _makeFilenameFromExpRun(self):
    if type(self.config.run) is not list:
      run = [self.config.run]
    else:
      run = self.config.run
    filenames = []
    for trun in run:
      tpath = '%s/hdf5'  %(self.config.experiment)
      tfile = '%s-r%04d.h5'  %(self.config.experiment,trun)
      filenames.append(os.path.join(self.config.dataPath,tpath,tfile))
    return filenames

  def _getFileStrategy(self):
    if self.config.fileNamesH5 == [] and not self.config.fileNamesIxp==[]:
      self.config.filestrategy = ['ixp']
    if not self.config.fileNamesH5 == [] and not self.config.fileNamesIxp==[]:
      if self.config.readCachedData:
	print "found ixp file in cache, will try to use it (override with readCachedData keyword)"
	self.config.filestrategy = ['h5','ixp']
      else:
	print "found ixp file in cache, will not be used, but might get overridden!!"
	self.config.filestrategy = ['h5']

    elif self.config.fileNamesH5 == [] and not self.config.fileNamesIxp==[]:
      self.config.filestrategy = ['ixp']
    elif not self.config.fileNamesH5 == [] and self.config.fileNamesIxp==[]:
      self.config.filestrategy = ['h5']

############## DATA TYPES #########################

class memdata(object):
  def __init__(self,name=None,input=None,scan=None,grid=None):
    self._name = name
    self.scan = scan
    self.grid = grid
    if input==None:
      raise Exception("memdata can not be initialized as no datasources were defined.")
    if hasattr(input,'isevtObj') and input.isevtObj:
      self._evtObj = input
      self._rdAllData = None
    elif hasattr(input,'__call__'):
      self._rdAllData = input
      self._evtObj = None
    else:
      self._evtObj = None
      self._rdAllData = None
      self._data,self._time = initmemdataraw(input)

      lens = [len(td) for td in self._data]
      self._filter = unravelScanSteps(np.arange(np.sum(lens)),lens)
      self._Nsteps = len(lens)
    self._expand = False
    self.interp  = interp(self)
    
    if not self._evtObj==None:
      self._evtObj.evt_cb = self._evtCall
  
  def _getFilteredData(self):
    if not hasattr(self,'_data'):
      if not self._rdAllData==None:
	self._data,self._time = initmemdataraw(self._rdAllData())
	lens = [len(td) for td in self._data]
	self._filter = unravelScanSteps(np.arange(np.sum(lens)),lens)
	self._Nsteps = len(lens)

    dat,stepsz = ravelScanSteps(self._data)
    return [dat[tf] for tf in self._filter]
  data = property(_getFilteredData)
  
  def _getFilteredTime(self):
    if not hasattr(self,'_time'):
      if not self._rdAllData==None:
	self._data,self._time = initmemdataraw(self._rdAllData())
	lens = [len(td) for td in self._data]
	self._filter = unravelScanSteps(np.arange(np.sum(lens)),lens)
	self._Nsteps = len(lens)
    dat,stepsz = ravelScanSteps(self._time)
    return [dat[tf] if len(tf)>0 else [] for tf in self._filter]
  time = property(_getFilteredTime)

  def __repr__(self):
    return "memdata associated to %s" % str(self._name) + '\n' + self._get_ana_str() 

  def _get_ana_str(self):
    ostr = ''
    ad = self.R
    hrange = np.percentile(ad,[5,95])

    formnum = lambda(num): '{:<9}'.format('%0.4g' %(num))
    for n in range(len(self)):
      ostr+='Step %04d:'%n + tools.hist_asciicontrast(self[n],bins=40,range=hrange,disprange=False) +'\n'
    ho = tools.hist_ascii(ad,range=hrange,bins=40)
    ostr+=ho.horizontal()
    return ostr
  
  def _evtCall(self):
    dat,time,stepNo,fileNo = self._evtObj.rdEvent()
    self._addEvt(dat,time,stepNo,fileNo)
  def _addEvt(self,data,time,stepNo,fileNo):
    pass
  
  def __len__(self):
    return len(self._filter)
  
  def _get_lens(self):
    return np.asarray([len(td) for td in self._filter])
  
  lens = property(_get_lens)
  
  def _dataSource(self):
    if not self._data==None:
      return 'data'
    else:
      if not self._rdStride==None:
	return 'rdStride'
      elif not self_evtObj==None:
	return 'evtObj'
      else:
	return None

  def ravel(self):
    return np.hstack(self.data)

  def filter(self,lims=None,inplace=False,perc=False):
    dat,stsz = ravelScanSteps(self.data)
    lims,filt = filter(dat,lims,perc=perc)
    filt = unravelIndexScanSteps(filt.nonzero()[0],stsz)
    if inplace:
      self._filter = filt
      return self
    else:
      tim,dum = ravelScanSteps(self.time)
      return memdata(input=[[dat[tf] for tf in filt],
                            [tim[tf] for tf in filt]],scan=self.scan,grid = self.grid)
  
  def digitize(self,bins=None,inplace = False):
    dat,stsz = ravelScanSteps(self.data)
    inds,bins = digitize(dat,bins)
    binrange = range(1,len(bins))

    filt = [(inds==tbin).nonzero()[0] for tbin in binrange]
    if inplace:
      self._filter = filt
      return self
    else:
      tim,dum = ravelScanSteps(self.time)
      scan = tools.dropObject(name='scan')
      if self._name==None:
	nname = ''
      else:
	nname = self._name+'_'
      scan[nname+'binedges'] = bins
      scan[nname+'bincenters'] = tools.histVecCenter(bins)

      return memdata(input=[[dat[tf] for tf in filt],
                            [tim[tf] for tf in filt]],scan=scan)


  def digitizeN(self,*args):
    return digitizeN(*args,target=self)
      #def digitize_new(self,*args,inplace=False):

  #def digitize_new(self,*args):
    #print args
    #dat,stsz = ravelScanSteps(self.data)
    #tim,stsz = ravelScanSteps(self.time)
    #dimension = len(args)
    #res_bins_dat = [dat]
    #res_bins_tim = [tim]
    #res_binvec = []
    #for tbins in args:
      #raise NotImplementeError
      #if isinstance(tbins,memdata):
	#print "Found memdata"
	#res_bins_dat_tmp = []
	#res_bins_tim_tmp = []
	#for bins_dat,bins_tim in zip(res_bins_dat,res_bins_tim):
	  #dum,tind = filterTimestamps(tbins.time,bins_tim)
	  #res_bins_tim_tmp.extend([bins_tim[ttind] for ttind in tind])
	  #res_bins_dat_tmp.extend([bins_dat[ttind] for ttind in tind])
	#res_bins_dat = res_bins_dat_tmp
        #res_bins_tim = res_bins_tim_tmp
	#res_binvec.append(tbin.scan[tbin.scan.keys()[0]])
      #else:
	#res_bins_dat_tmp = []
	#res_bins_tim_tmp = []
	#for bins_dat,bins_tim in zip(res_bins_dat,res_bins_tim):
	  #inds,bins = digitize(bins_dat,tbins)
	  ## continue here ToDo
	  #binrange = range(1,len(bins))
	  #tind = [(inds==tbin).nonzero()[0] for tbin in binrange]
	  #res_bins_tim_tmp.extend([bins_tim[ttind] for ttind in tind])
	  #res_bins_dat_tmp.extend([bins_dat[ttind] for ttind in tind])
	#res_bins_dat = res_bins_dat_tmp
        #res_bins_tim = res_bins_tim_tmp
      #return memdata(input=[res_bins_dat,res_bins_tim])


	#if inplace:
	  #self._filter = filt
	  #return self
	#else:
	  #tim,dum = ravelScanSteps(self.time)
	  #scan = tools.dropObject(name='scan')
	  #scan['binedges'] = bins
	  #scan['bincenters'] = tools.histVecCenter(bins)

	  #return memdata(input=[[dat[tf] for tf in filt],
				#[tim[tf] for tf in filt]],scan=scan)
  #def digitize_even_newer(self,*args):
    #dat,stsz = ravelScanSteps(self.data)
    #inds,bins = digitize(dat,bins)
    #binrange = range(1,len(bins))

    #filt = [(inds==tbin).nonzero()[0] for tbin in binrange]
    #if inplace:
      #self._filter = filt
      #return self
    #else:
      #tim,dum = ravelScanSteps(self.time)
      #scan = tools.dropObject(name='scan')
      #scan['binedges'] = bins
      #scan['bincenters'] = tools.histVecCenter(bins)

      #return memdata(input=[[dat[tf] for tf in filt],
                            #[tim[tf] for tf in filt]],scan=scan)
  def interpolate(self,other,type='mean',Nclosest=5):
    timenew = other.time
    o = getClosestEvents(self.time,timenew,Nclosest)
    sr,stepsz = ravelScanSteps(self.data)
    if type=='mean':
      op = np.asarray([np.mean(sr[tools.smartIdx(tsel)]) for tsel in o[0]])
      #de=bug
    odat = unravelScanSteps(op,o[1])
    return memdata(input=[odat,timenew],scan=other.scan)
  #return data(time=timenew,input=procObj,scan=self.scan)

  def ones(self):
    odat = [np.ones(len(dat)) for dat in self._data]
    return memdata(input=[odat,self.time],scan=self.scan,grid=self.grid)
  
  def _stepStatFunc(self,func,*args,**kwargs):
    data = np.asarray([func(d,*args,**kwargs) for d in self])
    if not self.grid==None:
      data = data.reshape(self.grid.shape)
    return data

  def mean(self,weights=None):
    if not weights==None:
      weights = self.ones()*weights
      selfdat = weights.ones()*self
      data =  np.asarray([np.average(d,weights=w) if len(w)>0 else np.nan for d,w in zip(selfdat,weights)])
      if not self.grid==None:
	data = data.reshape(self.grid.shape)
      return data
    else:
      return self._stepStatFunc(np.mean)
      
  def std(self,weights=None):
    if not weights==None:
      weights = self.ones()*weights
      selfdat = weights.ones()*self
      return np.asarray([tools.weighted_avg_and_std(d,w)[1] for d,w in zip(selfdat,weights)])
    else:
      return self._stepStatFunc(np.std)
  def median(self,weights=None):
    if not weights==None:
      weights = self.ones()*weights
      selfdat = weights.ones()*self
      return np.asarray([tools.weighted_median(d,w) for d,w in zip(selfdat,weights)])
    else:
      return self._stepStatFunc(np.median)
  def mad(self,weights=None):
    if not weights==None:
      weights = self.ones()*weights
      selfdat = weights.ones()*self
      return np.asarray([tools.weighted_mad(d,w) for d,w in zip(selfdat,weights)])
    else:
      return self._stepStatFunc(tools.mad)
  def sum(self):
    return self._stepStatFunc(np.sum)
  def count(self):
    return self._stepStatFunc(len)
  #def sqrt(self):
    #return np.asarray([np.sqrt(d) for d in self._data])
  def plot(self,weights=None):
    i = self.median(weights=weights).ravel()
    if not self.grid==None:
      x,y,im = self.grid.format(i)
      tools.imagesc(x,y,im.T)

    else:
      e = self.mad(weights=weights)/np.sqrt(self.lens)
      
      scanfields = self.scan.__dict__.keys()
      scanfields = [tf for tf in scanfields if not tf[0]=='_']
      x = self.scan.__dict__[scanfields[0]]
      tools.nfigure('memdata plot')
      pl.errorbar(x,i,yerr=e,fmt='.-')
      pl.xlabel(scanfields[0])


  def corrFilter(self,other,ratio=False):
    filt,lims = corrFilt(self,other,ratio=ratio)
    return self.filter(lims[0])

  R = property(ravel)


  def __getitem__(self,x):
    return self.data[x]
  
  # functions to make it feel like a datatype
  #def __add__(self,other):
    #return applyMemdataOperator(operator.add,self,other)
  #def __radd__(self,other):
    #return applyMemdataOperator(operator.add,self,other)
  #def __mul__(self,other):
    #return applyMemdataOperator(operator.mul,self,other)
  #def __rmul__(self,other):
    #return applyMemdataOperator(operator.mul,self,other)
  #def __div__(self,other):
    #return applyMemdataOperator(operator.div,self,other)
  #def __rdiv__(self,other):
    #return applyMemdataOperator(operator.div,self,other,isreverse=True)
  #def __truediv__(self,other):
    #return applyMemdataOperator(operator.truediv,self,other)
  #def __rtruediv__(self,other):
    #return applyMemdataOperator(operator.truediv,self,other,isreverse=True)
  #def __floordiv__(self,other):
    #return applyMemdataOperator(operator.floordiv,self,other)
  #def __rfloordiv__(self,other):
    #return applyMemdataOperator(operator.floordiv,self,other,isreverse=True)
  #def __mod__(self,other):
    #return applyMemdataOperator(operator.mod,self,other)
  #def __rmod__(self,other):
    #return applyMemdataOperator(operator.mod,self,other,isreverse=True)
  #def __sub__(self,other):
    #return applyMemdataOperator(operator.sub,self,other)
  #def __rsub__(self,other):
    #return applyMemdataOperator(operator.sub,self,other,isreverse=True)
  #def __pow__(self,other):
    #return applyMemdataOperator(operator.pow,self,other)
  #def __rpow__(self,other):
    #return applyMemdataOperator(operator.pow,self,other,isreverse=True)

  #def __and__(self,other):
    #return applyMemdataOperator(operator.and_,self,other)
  #def __rand__(self,other):
    #return applyMemdataOperator(operator.and_,self,other)
  #def __or__(self,other):
    #return applyMemdataOperator(operator.or_,self,other)
  #def __ror__(self,other):
    #return applyMemdataOperator(operator.or_,self,other)
  #def __xor__(self,other):
    #return applyMemdataOperator(operator.xor,self,other)
  #def __rxor__(self,other):
    #return applyMemdataOperator(operator.xor,self,other)
  #def __le__(self,other):
    #return applyMemdataOperator(operator.le,self,other)
  #def __lt__(self,other):
    #return applyMemdataOperator(operator.lt,self,other)
  #def __eq__(self,other):
    #return applyMemdataOperator(operator.eq,self,other)
  #def __ne__(self,other):
    #return applyMemdataOperator(operator.ne,self,other)
  #def __ge__(self,other):
    #return applyMemdataOperator(operator.ge,self,other)
  #def __gt__(self,other):
    #return applyMemdataOperator(operator.gt,self,other)
  # try if works
  def __add__(self,other):
    return applyDataOperator(operator.add,self,other)
  def __radd__(self,other):
    return applyDataOperator(operator.add,self,other)
  def __mul__(self,other):
    return applyDataOperator(operator.mul,self,other)
  def __rmul__(self,other):
    return applyDataOperator(operator.mul,self,other)
  def __div__(self,other):
    return applyDataOperator(operator.div,self,other)
  def __rdiv__(self,other):
    return applyDataOperator(operator.div,self,other,isreverse=True)
  def __truediv__(self,other):
    return applyDataOperator(operator.truediv,self,other)
  def __rtruediv__(self,other):
    return applyDataOperator(operator.truediv,self,other,isreverse=True)
  def __floordiv__(self,other):
    return applyDataOperator(operator.floordiv,self,other)
  def __rfloordiv__(self,other):
    return applyDataOperator(operator.floordiv,self,other,isreverse=True)
  def __mod__(self,other):
    return applyDataOperator(operator.mod,self,other)
  def __rmod__(self,other):
    return applyDataOperator(operator.mod,self,other,isreverse=True)
  def __sub__(self,other):
    return applyDataOperator(operator.sub,self,other)
  def __rsub__(self,other):
    return applyDataOperator(operator.sub,self,other,isreverse=True)
  def __pow__(self,other):
    return applyDataOperator(operator.pow,self,other)
  def __rpow__(self,other):
    return applyDataOperator(operator.pow,self,other,isreverse=True)

  def __and__(self,other):
    return applyDataOperator(operator.and_,self,other)
  def __rand__(self,other):
    return applyDataOperator(operator.and_,self,other)
  def __or__(self,other):
    return applyDataOperator(operator.or_,self,other)
  def __ror__(self,other):
    return applyDataOperator(operator.or_,self,other,isreverse=True)
  def __xor__(self,other):
    return applyDataOperator(operator.xor,self,other)
  def __rxor__(self,other):
    return applyDataOperator(operator.xor,self,other,isreverse=True)
  def __le__(self,other):
    return applyDataOperator(operator.le,self,other)
  def __lt__(self,other):
    return applyDataOperator(operator.lt,self,other)
  def __eq__(self,other):
    if other==None:
      return False
    return applyDataOperator(operator.eq,self,other)
  def __ne__(self,other):
    if other==None:
      return True
    return applyDataOperator(operator.ne,self,other)
  def __ge__(self,other):
    return applyDataOperator(operator.ge,self,other)
  def __gt__(self,other):
    return applyDataOperator(operator.gt,self,other)
  
  def __invert__(self):
    return applyDataOperator(operator.invert,self)

class data(object):
  def __init__(self,name=None,time=None,input=None,scan=None,ixpAddressData=None,parent=None,grid=None):
    self.name = name
    self.scan = scan
    self.grid = grid
    self._ixpAddressData = ixpAddressData
    self._parent = parent
    #if input==None:
      #raise Exception("data can not be initialized as no datasources were defined.")
    if hasattr(input,'isevtObj') and input.isevtObj:
      self._procObj = None
      self._evtObj = input
      self._rdStride = None
    elif hasattr(input,'__call__'):
      self._procObj = None
      self._rdStride = input
      self._evtObj = None
    elif type(input) is dict:
      self._procObj = input
      self._evtObj = None
      self._rdStride = self._rdFromObj
    elif type(input) is h5py.highlevel.Group:
      self._procObj = None
      self._evtObj = None
      self._rdStride = self._rdFromIxp
    else:
      self._procObj = None
      self._evtObj = None
      self._rdStride = None
    if not self._ixpAddressData is None:
      self._rdStride = self._rdFromIxp

    if time==None:
      self._time = []
      print "Warning: Missing timestamps for data instance,\n--> data needs to be appended."
    else:
      if hasattr(time,'isevtObj') and input.isevtObj:
	self._evtObjtime = input
	self._rdTime = None
      elif hasattr(time,'__call__'):
	self._rdTime = input
	self._evtObj = None
      else:
	self._time = time
	self._evtObj = None
	self._rdTime = None
    lens = [len(td) for td in self._time]
    self._filter = unravelScanSteps(np.arange(np.sum(lens)),lens)
    self._Nsteps = len(lens)
    self._lens = lens
    self.evaluate = Evaluate(self)
    self._cacheLast = False                                                          
    self._lastCache = None
    self._showMe = False
    self._showMeName = None

  def _showme(self,name=None):
    self._showMe = True
    self._showMeName = name


  def _deleteIxp(self,force=False):
    ixp,path = self._getIxpHandle()
    if path in ixp.fileHandle and not force:
      deleteit = 'y' == raw_input("Dataset %s exists in ixp file, would you like   to delete it? (y/n) "%path)
      if deleteit:
	del ixp.fileHandle[path]
      else:
	return

  def _appendData(self,input_data):
    assert np.iterable(input_data), "input data has to be iterable"
    assert (type(input_data) is list) or (type(input_data) is tuple)  , "input_data has to be either list or tuple"
    data = input_data[0]
    time = input_data[1]
    if len(input_data)>2:
      steps = input_data[2]
    else:
      steps = None
    ixp,path = self._getIxpHandle()
    grp = ixp.fileHandle.require_group(path)
    if (type(data) is list) and (type(time) is list):
      self._time.extend(time)
      for tdata,ttime in zip(data,time):
	ixp.appendData(grp,(tdata,ttime))

    self._ixpAddressData = grp['data']
    self._rdStride = self._rdFromIxp
    lens = [len(td) for td in self._time]
    self._filter = unravelScanSteps(np.arange(np.sum(lens)),lens)
    self._Nsteps = len(lens)
    self._lens = lens
  
  def __len__(self):
    return len(self._lens)

  def lens(self):
    return np.asarray([len(td) for td in self._filter])
  
  def ones(self):
    odat = [np.ones(len(tim)) for tim in self.time]
    return memdata(input=[odat,self.time],scan=self.scan)

  def __repr__(self):
      return "`data` object %s, %d steps, %d events per step" % (self.name, self.__len__(),np.median(self._lens))

  def _getFilteredTime(self):
    if not hasattr(self,'_time'):
      if not self._rdAllData==None:
	self._data,self._time = initmemdataraw(self._rdAllData())
	lens = [len(td) for td in self._data]
	self._filter = unravelScanSteps(np.arange(np.sum(lens)),lens)
	self._Nsteps = len(lens)
    dat,stepsz = ravelScanSteps(self._time)
    return [dat[tf] for tf in self._filter]
  time = property(_getFilteredTime)
  ## data object (self)manipulation
  #def _getFromSelf(self,what):
    #""" get data from object for example: d.ipm2.get("_file0.step3.channel") """
    #return tools.getFromObj(self,what)
  #def _existsInSelf(self,what):
    #return tools.existsInObj(self,what)
  #def _addToSelf(self,what,value,overWrite=True):
    #return tools.addToObj(self,what,value,overWrite=overWrite)

  def _memIterate(self,steps=slice(None),evts=slice(None),memFrac=0.1):
    steps = tools.itemgetToIndices(steps,len(self._filter))
    evts = [tools.itemgetToIndices(evts,len(self._filter[step])) for step in steps]
    sizeEvt = None
    if not hasattr(self,'_mem'):
      self._mem = mem()
    indchunks = []
    for stepNo,step in enumerate(steps):
      tevts = evts[stepNo]

      tlen  = len(tevts)
      nextEvtNo = 0
      stepChunks = []
      if not hasattr(self,'_sizeEvt') or (self._sizeEvt is None):
	usedmem=0
	thisstuff = self
	done = False
	gotshape = False
	while not done:
          tdat = thisstuff[step,tevts[0]]
	  if not gotshape:
	    thisshape = np.shape(tdat[0][0])
	    gotshape = True
	  usedmem += tdat[0].nbytes
	  print usedmem
	  if thisstuff._procObj==None:
	    done = True
	  else:
	    args = thisstuff._procObj['args']
	    isdatainst = np.asarray([isinstance(targ,data) for targ in args])
	    if np.sum(isdatainst)>0:
	      thisstuff = args[isdatainst.nonzero()[0][0]]
	    elif 'ixppyInput' in thisstuff._procObj.keys():
	      ixpIP = thisstuff._procObj['ixppyInput'][0][0]
	      thisstuff = ixpIP


        self._sizeEvt = dict(bytes=usedmem , shape=thisshape)
      Nread = np.floor(self._mem.updatefree()/8/self._sizeEvt['bytes']*memFrac)
      while nextEvtNo<tlen:
        stepChunks.append(tevts[nextEvtNo:min(tlen-1,nextEvtNo+Nread)])
        nextEvtNo+=Nread
      indchunks.append((stepNo,stepChunks))
    return indchunks
  
  def chunks(self,steps=slice(None),evts=slice(None),memFrac=0.1):
    indchunks = self._memIterate(steps,evts,memFrac)
    out = []
    stepslen = len(indchunks)
    for sNo,step in enumerate(indchunks):
      stepd = []
      for chunk in step:
        td = pycopy.copy(self)
        td._filter = pycopy.deepcopy(td._filter)
        #td._filter = [[],] * stepslen
        #td._filter[sNo] = np.sum(td._lens[:sNo])+chunk
        td._filter = [np.sum(td._lens[:sNo])+chunk]
        stepd.append(td)
      out.append(stepd)
    return out

  def _getIxpHandle(self):
    ixpSeed = None
    path = ''
    present = self
    while ixpSeed is None:
      tp = present._parent
      name = [name for name,tid in tp.__dict__.iteritems() if id(tid)==id(present)][0]
      if hasattr(tp,'_ixpHandle'):
	ixpSeed = tp._ixpHandle
	path = '/dataset/'+name+'/'+path 
      else:
	present = tp
	path = name+'/'+path
    return ixpSeed,path

  def get_memdata(self):
    data = self[:,:]
    nel = np.shape(data[0])[1]
    memdat = []
    for n in range(nel):
      tmemdat = []
      for step in data:
	tmemdat.append(step[:,n])
      memdat.append(memdata(input=[tmemdat,self.time],scan=self.scan,grid=self.grid))
    return memdat

  #def _copy(self):
  ##data(time=self.time,
  ##def __init__(self,name=None,time=None,input=None,scan=None,ixpAddress=None,parent=None):
  ##self.name = name
  ##self.scan = scan
  ##self._ixpAddress = ixpAddress
  ##self._parent = parent
    #if self._proc
    #if input==None:
      #raise Exception("data can not be initialized as no datasources were defined.")
    #if hasattr(input,'isevtObj') and input.isevtObj:
      #self._procObj = None
      #self._evtObj = input
      #self._rdStride = None
    #elif hasattr(input,'__call__'):
      #self._procObj = None
      #self._rdStride = input
      #self._evtObj = None
    #elif type(input) is dict:
      #self._procObj = input
      #self._evtObj = None
      #self._rdStride = self._rdFromObj
    #else:
      #self._procObj = None
      #self._evtObj = None
      #self._rdStride = None

    #if time==None:
      #print Exception("Warning: Missing timestamps for data instance!")
    #if hasattr(time,'isevtObj') and input.isevtObj:
      #self._evtObjtime = input
      #self._rdTime = None
    #elif hasattr(time,'__call__'):
      #self._rdTime = input
      #self._evtObj = None
    #else:
      #self._time = time
      #self._evtObj = None
      #self._rdTime = None
    #lens = [len(td) for td in self._time]
    #self._filter = unravelScanSteps(np.arange(np.sum(lens)),lens)
    #self._Nsteps = len(lens)
    #self._lens = lens
 
    
  def interpolate(self,timenew,type='mean',Nclosest=5):
    o = getClosestEvents(self.time,timenew,Nclosest)
    if type=='mean':
      def func(crossInput=[]):
        return np.mean(crossInput[0],axis=0)

    procObj = dict(
	func=func,
	ixppyInput=[(self,o[0])],
	isPerEvt=False,
	args=[],
	kwargs=dict(),
	isCrossEvent=True,
	nargSelf=0)
    return data(time=timenew,input=procObj,scan=self.scan)

  def __getitem__(self,x):
    n = self._Nsteps
    lens = [len(tfilt) for tfilt in self._filter]
    x = tools.iterfy(x)
    if len(x)==1:
      if n==1:
        stepInd = [0]
	#if len(x[1])==0:
	  #return []
        evtInd = [tools.itemgetToIndices(x[0],lens[tind],boolean=True) for tind in stepInd]
      else:
	raise IndexError('More than one scanstep [scanstep,event] index pair required!')
    elif len(x)==2:
      stepInd = tools.itemgetToIndices(x[0],n)
      #if len(x[1])==0:
	#return []
      evtInd = [tools.itemgetToIndices(x[1],lens[tind]) for tind in stepInd]
    #print stepInd,evtInd
    return self._getStepsShots(stepInd,evtInd)
    

  def _getStepsShots(self,stepInd,evtInd):
    if self._showMe:
      t0 = time.time()
      print "%s is requested for stepInd %s and evtInd %s"%(self._showMeName,stepInd,evtInd)
    if self._cacheLast:
      if self._lastCache is None:
        readit = True
      else:
        if stepInd==self._lastCache['stepInd']\
            and np.asarray([len(teI)==len(tceI) for teI,tceI in zip(evtInd,self._lastCache['evtInd'])]).all()\
	    and np.asarray([(teI==tceI).all() for teI,tceI in zip(evtInd,self._lastCache['evtInd'])]).all():
	  dat = self._lastCache['dat']
	  readit=False
	  print "recycle"
	else:
	  readit=True
    else:
      readit=True

    if readit:
      lens = [len(tfilt) for tfilt in self._filter]
      inds = ravelIndexScanSteps(evtInd,lens,stepNo = stepInd)
      indsdat = np.hstack(self._filter).argsort()[inds]
      indsdat_sorted = np.sort(indsdat)
      indsdat_read = unravelIndexScanSteps(indsdat_sorted,self._lens)
      stepInd_read = [n for n in range(len(indsdat_read)) if len(indsdat_read[n])>0]
      evtInd_read  = [indsdat_read[n] - int(np.sum(self._lens[:n])) for n in range(len(indsdat_read)) if len(indsdat_read[n])>0]

      dat = [self._rdStride(step,tevtInd)[0] for step,tevtInd in zip(stepInd_read,evtInd_read)]
      if not dat==[]:
	dat = np.concatenate(dat)
	ind_resort = unravelIndexScanSteps(inds[inds.argsort()],lens,replacements=inds.argsort())
	ind_resort = [tools.smartIdx(l) for l in ind_resort if len(l)>0]
	#print [tools.smartIdx(l) for l in ind_resort if len(l)>0]
	dat = [dat[tind_resort,...] for tind_resort in ind_resort]
	if 0 in lens:
	  zinds = (np.asarray(lens)==0).nonzero()[0]
	  for zind in zinds:
	    if zind in stepInd:
	      dat.insert(zind,[])

      self._lastCache = dict(dat=dat,stepInd=stepInd,evtInd=evtInd)	

    if self._showMe:
      print "...this took %4g seconds." %(time.time()-t0)
    return dat

    #if len(stepInd)>1:
    #return [self._rdStride(step,tevtInd) for step,tevtInd in zip(stepInd,evtInd)]
  def _applyfun(self,procObj,stepStride,eventStride):
    #TODO: 
    pass
  def _rdFromObj(self,step, evtInd):
    if self._procObj.has_key('isCrossEvent') and self._procObj['isCrossEvent']:
      ret = applyCrossFunction(self._procObj['func'],
	                       ixppyInput=self._procObj['ixppyInput'],
                               time=self.time,
			       args=self._procObj['args'],
			       kwargs=self._procObj['kwargs'],
			       stride = [step,evtInd])

    else:
      ret = applyFunction(self._procObj['func'],
			   self._procObj['args'],
			   self._procObj['kwargs'],
			   stride = [step,evtInd],
			   isPerEvt = self._procObj['isPerEvt'],
			   InputDependentOutput=True,
			   NdataOut=1,NmemdataOut=0, picky=False)
      if type(ret) is not tuple:
	ret = (ret,)

    #return [tret[self._procObj['nargSelf']] for tret in ret]
    return ret[self._procObj['nargSelf']]

  def _rdFromIxp(self,step,evtInd):
    if len(evtInd)==1:
      return [np.asarray([(self._ixpAddressData['#%06d'%step][np.asarray(evtInd),...])])]
    else:
      return [np.atleast_2d(self._ixpAddressData['#%06d'%step][np.asarray(evtInd),...])]
  def _isIxp(self):
    return self._rdStride == self._rdFromIxp

  def __add__(self,other):
    return applyDataOperator(operator.add,self,other)
  def __radd__(self,other):
    return applyDataOperator(operator.add,self,other)
  def __mul__(self,other):
    return applyDataOperator(operator.mul,self,other)
  def __rmul__(self,other):
    return applyDataOperator(operator.mul,self,other)
  def __div__(self,other):
    return applyDataOperator(operator.div,self,other)
  def __rdiv__(self,other):
    return applyDataOperator(operator.div,self,other,isreverse=True)
  def __truediv__(self,other):
    return applyDataOperator(operator.truediv,self,other)
  def __rtruediv__(self,other):
    return applyDataOperator(operator.truediv,self,other,isreverse=True)
  def __floordiv__(self,other):
    return applyDataOperator(operator.floordiv,self,other)
  def __rfloordiv__(self,other):
    return applyDataOperator(operator.floordiv,self,other,isreverse=True)
  def __mod__(self,other):
    return applyDataOperator(operator.mod,self,other)
  def __rmod__(self,other):
    return applyDataOperator(operator.mod,self,other,isreverse=True)
  def __sub__(self,other):
    return applyDataOperator(operator.sub,self,other)
  def __rsub__(self,other):
    return applyDataOperator(operator.sub,self,other,isreverse=True)
  def __pow__(self,other):
    return applyDataOperator(operator.pow,self,other)
  def __rpow__(self,other):
    return applyDataOperator(operator.pow,self,other,isreverse=True)

  def __and__(self,other):
    return applyDataOperator(operator.and_,self,other)
  def __rand__(self,other):
    return applyDataOperator(operator.and_,self,other)
  def __or__(self,other):
    return applyDataOperator(operator.or_,self,other)
  def __ror__(self,other):
    return applyDataOperator(operator.or_,self,other)
  def __xor__(self,other):
    return applyDataOperator(operator.xor,self,other)
  def __rxor__(self,other):
    return applyDataOperator(operator.xor,self,other)
  def __le__(self,other):
    return applyDataOperator(operator.le,self,other)
  def __lt__(self,other):
    return applyDataOperator(operator.lt,self,other)
  def __eq__(self,other):
    if other==None:
      return False
    return applyDataOperator(operator.eq,self,other)
  def __ne__(self,other):
    if other==None:
      return True
    return applyDataOperator(operator.eq,self,other)
  def __ge__(self,other):
    return applyDataOperator(operator.ge,self,other)
  def __gt__(self,other):
    return applyDataOperator(operator.gt,self,other)

  def mean(self):
    allchunks = self._memIterate()
    out = []
    # progress bar
    Nevtot = np.sum([ np.sum([len(tchunk) for tchunk in step]) for stepNo,step in allchunks])
    processedevents = 0
    widgets = ['Evaluating mean: ', pb.Percentage(), ' ', pb.Bar(),' ', pb.ETA(),'  ']
    pbar = pb.ProgressBar(widgets=widgets, maxval=Nevtot).start()
    for stepNo,step in allchunks:
      totlen = np.sum([len(x) for x in step])
      tout = 0
      for chunk in step:
	if len(chunk)>0:
          tout += np.nansum(self[stepNo,np.ix_(chunk)][0],axis=0)
        processedevents += len(chunk)
	pbar.update(processedevents)
      out.append(tout/totlen)
    pbar.finish()
    return out


Data = data
Memdata = memdata

########## DATA TYPE OPERATION ###########################

def applyOperator(optr,a,b,isreverse=False):
  a = tools.iterfy(a)
  b = tools.iterfy(b)

  res = []
  if not isreverse:
    if len(a)==len(b):
      for ta,tb in zip(a,b):
          res.append(optr(ta,tb))
    else:
      for ta in a:
        res.append(optr(ta,b))
  else:
    if len(a)==len(b):
      for ta,tb in zip(a,b):
          res.append(optr(tb,ta))
    else:
      for ta in a:
        res.append(optr(b,ta))
  return res

def expandMemdata(a,b):
  aex = a._expand
  bex = b._expand
  if aex and bex:
    raise Exception("Can not expand both ixppy.memdata instances in one operation!")
  elif not aex and not bex:
    return a,b
  else:
    if aex:
      pass
      
def applyMemdataOperator(optr,a,b,isreverse=False):
  amem = isinstance(a,memdata) 
  bmem = isinstance(b,memdata) 

  if amem and bmem:
    a,b = expandMemdata(a,b)
    #a,b = interpMemdata(a,b)
    #a,b = interpStepMemdata(a,b)

    ai,bi = filterTimestamps(a.time,b.time)
    # TODO: check if ressource expensive
    ar = np.hstack(a.data)
    br = np.hstack(b.data)
    ri,rstepsizes = ravelScanSteps(ai)
    resdat = optr(ar[ri],br[np.hstack(bi)])
    restim = np.hstack(a.time)[ri]
    resdat = unravelScanSteps(resdat,rstepsizes)
    restim = unravelScanSteps(restim,rstepsizes)
    scan   = a.scan 
  elif amem or bmem:
    if amem:
      adat = a.data
      restim = a.time
      bdat = b
      scan   = a.scan 
    elif bmem:
      bdat = b.data
      restim = b.time
      adat = a
      scan   = a.scan 
    adat = tools.iterfy(adat)
    bdat = tools.iterfy(bdat)
    resdat = []
    if not isreverse:
      if len(adat)==len(bdat):
	for ta,tb in zip(adat,bdat):
	    resdat.append(optr(ta,tb))
      else:
	for ta in adat:
	  resdat.append(optr(ta,bdat))
    else:
      if len(adat)==len(bdat):
	for ta,tb in zip(adat,bdat):
	    resdat.append(optr(tb,ta))
      else:
	for ta in adat:
	  resdat.append(optr(bdat,ta))
  
  return memdata(input=[resdat,restim],scan=scan)
      
def applyDataOperator(optr,a,b=None,isreverse=False):
  if b==None:
    args=[a]
  else:
    if not isreverse:
      args = [a,b]
    else:
      args = [b,a]
  return applyFunction(optr,args,dict(),InputDependentOutput=True, NdataOut=1,NmemdataOut=0, picky=False, isPerEvt=False, outputtypes=None)

def _applyFun(func,a):
  res = ixppyList()
  for ta in a:
    res.append(func(ta))
  return res

def applyFunction(func,ipargs,ipkwargs,InputDependentOutput=True, KWignore=None, NdataOut=0,NmemdataOut=0, picky=False, isPerEvt=False, stride=None, outputtypes=None, forceCalculation=False):
  """ rules: 
  - if data output, no other output possible, as output is not calculated. 
  - all event dependent arguments have to be passed as memdata or data instances
  - first arg has timestamp sort priority, kwargs after args, kwargs are inter-
    preted in "sort- order". This priority order also holds for scanvec data.
  - When InputDependentOutput is true (default) the output expects:
    a) one data instance if any data insctance in input
    b) one memdata instance if memdata instance in input and no data instance
    otherwise NmemdataOut and NdataOut need to be specified. All other output 
    comes after ixppy instances if avaiable.
  """
  if outputtypes==None:
    outputtypes = NdataOut*['data'] + NmemdataOut*['memdata']
  ##### Filter timestampsin order, args first, kwargs after sorted keys
  allobjects = [(arg,0,argno) for argno,arg in enumerate(ipargs) if (isinstance(arg,data) or isinstance(arg,memdata))]
  kwkeys = ipkwargs.keys()
  kwkeys.sort()
  if KWignore==None:
    KWignore = []
  else:
    KWignore = iterfy(KWignore)
  for keyno,key in enumerate(kwkeys):
    if not key in KWignore:
      if (isinstance(ipkwargs[key],data) or isinstance(ipkwargs[key],memdata)):
	allobjects.append((ipkwargs[key],1,keyno))
  if not allobjects==[]:
    scan = allobjects[0][0].scan
    grid = allobjects[0][0].grid
    
  rtimes = get_common_timestamps([ao[0] for ao in allobjects])
  # get also other arguments in seperate list for later use in length analysis if needed
  if not picky and not rtimes==None:
    # get other input
    otherargs = [(arg,0,argno) for argno,arg in enumerate(ipargs) if ( not isinstance(arg,data) and not isinstance(arg,memdata))]
    for keyno,key in enumerate(kwkeys):
      if ( not isinstance(arg,data) and not isinstance(arg,memdata)):
	otherargs.append((ipkwargs[key],1,keyno))
    otherlens = [(len(toa),iskey,argno) for toa,iskey,argno in otherargs if (type(toa) is not str and np.iterable(toa))]
    if not otherlens==[]:
      lens,lensiskey,lensargno = zip(*otherlens)
      #tequallengths = np.logical_and(np.asarray(lens)==len(rtimes),np.logical_not(np.asarray(lens)==1))
      tequallengths = np.logical_and(np.asarray(lens)==len(rtimes),True)
      otherip = [otherargs[eqi] for eqi in tequallengths.nonzero()[0]]
    else:
      otherip = []
  else:
    otherip = []

  ############ generate output ############
  # case of data instance in input, at the moment seems like data might come out, but this has to be thought about more.
  if np.asarray([isinstance(to[0],data) for to in allobjects]).any() or forceCalculation:
    if 'data' in outputtypes and stride==None:
      # this is the normal case to make an object that will act upon call
      output = []
      for nargSelf in (np.asarray(outputtypes)=='data').nonzero()[0]:
        procObj = dict(func=func,args=ipargs,kwargs=ipkwargs,nargSelf=nargSelf,isPerEvt=isPerEvt)
        output.append(data(time=rtimes,input=procObj,scan=scan))
      if len(output)>1:
	output = tuple(output)
      else:
	output = output[0]
    # case mainly when executes as data procObj
    elif ('data' in outputtypes and not stride==None) or forceCalculation:
      # in force case and when stride is not given
      if (stride == None) and forceCalculation:
	# find smallest chunk size, will go for that...
	dataIPchunkings = [o._memIterate() for o,dum,dum in allobjects if isinstance(o,data)]
	chunksize = np.min([np.min([len(tchunk[1]) for tchunk in tchunks]) for tchunks in dataIPchunkings])
	# Get step stride and event strides
	stepstride = range(len(rtimes))
	eventstrides = []
	for trtimes in rtimes:
	  evlst = range(len(trtimes))
	  eventstrides.append([ evlst[i:i+chunksize] for i in range(0, len(evlst), chunksize) ])
      else:
	stepstride = tools.iterfy(stride[0])
	eventstrides = [stride[1]]*len(stepstride)

      ixppyip = []
      ixppytype = []
      for o,iskey,argind in allobjects: 
	ir,io      = filterTimestamps(rtimes,o.time)
	
	if isinstance(o,data):
	  ixppytype.append('data')
	  eventstrides_ravel = [np.ravel(tevs) for tevs in eventstrides]

	  io = getStepShotsFromIndsTime(io,o.time,stride=eventstrides_ravel,getIncludedSteps=True)
	if isinstance(o,memdata):
	  ixppytype.append('memdata')
	ixppyip.append(([o,io,iskey,argind]))

      # generate input structures
      output_list = []

      for stepstrideNo,step in enumerate(stepstride):
	eventstride = eventstrides[stepstrideNo]
	if type(eventstride[0]) is list:
	  print "this chunking doesn't work yet, taking first chunk only"
	  eventstride = eventstride[0]

	targs   = list(pycopy.copy(ipargs))
	tkwargs = pycopy.copy(ipkwargs)
	for o,io,k,i in ixppyip:
	  if not k:
	    if isinstance(o,memdata):
	      odat,stpsz = ravelScanSteps(o.data)
	      targs[i] = odat[io[step]][eventstride]
	    else:
	      
	      #tmp = o._getStepsShots(io[step][0],io[step][1])[0]
	      
	      io,inclsteps = io
	      tmp = np.concatenate(o._getStepsShots(io[(inclsteps==step).nonzero()[0]][0],io[(inclsteps==step).nonzero()[0]][1]),axis=0)
	      if not isPerEvt and ('memdata' in ixppytype):
	        trnspsorder = range(np.rank(tmp))
                trnspsorder = trnspsorder[1:]+[trnspsorder[0]]
                tmp         = tmp.transpose(trnspsorder)
              targs[i] = tmp
	      #raise NotImplementedError('Use the source, luke!')

	  else:
	    tkwargs[kwkeys[i]]  = o[step][eventstride]
	    if isinstance(o,memdata):
	      odat,stpsz = ravelScanSteps(o.data)
	      tkwargs[kwkeys[i]] = odat[io[step]][eventstride]
	    else:
              tmp = o._getStepsShots(io[step][0],io[step][1])[0]
	      if not isPerEvt and ('memdata' in ixppytype):
                trnspsorder = range(np.rank(tmp))
                trnspsorder = trnspsorder[1:]+[trnspsorder[0]]
                tmp = tmp.transpose(trnspsorder)
              tkwargs[kwkeys[i]] = tmp
	for o,k,i in otherip:
	  if not k:
	    targs[i] = o[step]
	  else:
	    tkwargs[kwkeys[i]]  = o[step]
	if isPerEvt:
	  # TODO: in case func works only for single shot
	  tret = []
	  for nevt in range(len(eventstride)):
            stargs = pycopy.copy(targs)
            stkwargs = pycopy.copy(tkwargs)
	    for o,io,k,i in ixppyip:
	      if not k:
                stargs[i] = targs[i][nevt]
	      else:
                stkwargs[kwkeys[i]] = tkwargs[kwkeys[i]][nevt]
            stret = func(*stargs,**stkwargs)
	    if type(stret) is not tuple: 
              stret = (stret,)
            tret.append(stret)
	  tret = tuple(zip(*tret))



	else:

	  tret = func(*targs,**tkwargs)
	  if not type(tret) is tuple:
            tret = (tret,)
	if not isPerEvt and ('memdata' in ixppytype):
	  tret = list(tret)
	  for ono,ttret in enumerate(tret):
	    rnk = np.rank(ttret)
	    if rnk>1:
	      trnspsorder = range(rnk)
	      trnspsorder = [trnspsorder[-1]]+trnspsorder[:-1]
	      tret[ono]   = ttret.transpose(trnspsorder)
	  tret = tuple(tret)
	output_list.append(tret)
    ############ interprete output automatically find memdata candidates ###########
      
      output_list = zip(*output_list)
      output_ismemdata = [opNo  for opNo,ao in enumerate(output_list) \
	  if [len(rtime) for rtime in rtimes] == [len(tools.iterfy(tao)) for tao in ao]]
      for n,top in enumerate(output_list):
	if n in output_ismemdata:
	  output_list[n] = memdata(input=[top,rtimes],scan=scan,grid=grid)
	#elif top.count(top[0])==len(top):
	  #output_list[n] = top[0]
	#elif len(top)>1:
        else:
	  output_list[n] = list(top)
      if len(output_list)>1:
	output = tuple(output_list)
      else:
	output = output_list[0]
      

  # "Easy" case where everything fits in memory, no data instance in input.
  elif np.asarray([isinstance(to[0],memdata) for to in allobjects]).any():
    ixppyip = []
    for o,iskey,argind in allobjects: #assuming here that data instances are out!
      ir,io      = filterTimestamps(rtimes,o.time)
      odat,stpsz = ravelScanSteps(o.data)
      ixppyip.append(([odat[tio] for tio in io],iskey,argind))
    # generate input structures
    output_list = []
    for step in range(len(rtimes)):
      targs   = list(pycopy.copy(ipargs))
      tkwargs = pycopy.copy(ipkwargs)
      for o,k,i in ixppyip + otherip:
	if not k:
	  targs[i] = o[step]
	else:
	  tkwargs[kwkeys[i]]  = o[step]
      if isPerEvt:
	# TODO: in case func works only for single shot
	pass
      else:
	tret = func(*targs,**tkwargs)
      if type(tret) is not tuple: 
	tret = (tret,)
      output_list.append(tret)
  ############ interprete output automatically find memdata candidates ###########
    output_list = zip(*output_list)
    # check for equal output and find candidates for data/memdata instances
    output_ismemdata = [opNo  for opNo,ao in enumerate(output_list) \
	if [len(rtime) for rtime in rtimes] == [len(tools.iterfy(tao)) for tao in ao]]
    for n,top in enumerate(output_list):
      if n in output_ismemdata:
	output_list[n] = memdata(input=[top,rtimes],scan=scan,grid=grid)
      elif top.count(top[0])==len(top):
	output_list[n] = top[0]
      elif len(top)>1:
	output_list[n] = list(top)
    if len(output_list)>1:
      output = tuple(output_list)
    else:
      output = output_list[0]
  else:
    #pass
    output = func(*ipargs,**ipkwargs)
  return output

def wrapFunc(func,InputDependentOutput=True, NdataOut=1,NmemdataOut=0, picky=False, isPerEvt=False, stride=None):
  """ rules: 
  - if data output, no other output possible, as output is not calculated. 
  - all event dependent arguments have to be passed as memdata or data instances
  - first arg has timestamp sort priority, kwargs after args, kwargs are inter-
    preted in "sort- order". This priority order also holds for scanvec data.
  - When InputDependentOutput is true (default) the output expects:
    a) one data instance if any data insctance in input
    b) one memdata instance if memdata instance in input and no data instance
    otherwise NmemdataOut and NdataOut need to be specified. All other output 
    comes after ixppy instances if avaiable.
  """

  @wraps(func,assigned=('__name__', '__doc__'))
  def wrapper(*args,**kwargs):
    return applyFunction(func,args,kwargs,InputDependentOutput=True, NdataOut=NdataOut,NmemdataOut=NmemdataOut, picky=picky, isPerEvt=isPerEvt, stride=None)
  return wrapper
    
def applyCrossFunction(func,ixppyInput=[], time=None, args=None,kwargs=None, stride=None, outputtypes=None):
  """ important keywords:
   time: list of timestamps which characterize the event
   ixppyInput: a list of (ixppytype,coordination) tuples 
               (ixppytype is memdata or data instance)
   other args andkwargs are possible. 
   
   func takes
   crossInput: list of elements to calculate stuff 

  """
  lens = [len(ttime) for ttime in time]
  inds = []
  for sNO,tlen in enumerate(lens):
    inds.append(np.int32(np.arange(tlen)+np.sum(lens[:sNO])))

  if stride is None:
    stride = [range(len(lens)),'all']

  output = []
  stepstride = tools.iterfy(stride[0])
  for step in stepstride:
    # get the indices from this ixppy type instance and step
    tstride = stride[1]
    if tstride=='all':
      tstride = range(lens[step])
    nind = [inds[step][tind] for tind in tstride]
    # Read necessary data
    package = []
    for tixpIp in ixppyInput:
      #unpacking...
      tdata = tixpIp[0]
      coo = tixpIp[1]
      # groups of inds per evt in this step
      xind = [coo[tnind] for tnind in nind]
      #translating this into steps and indices in tdata
      uxind,repack = np.unique(np.hstack(xind),return_inverse=True)
      tstride = getStepShotsFromIndsTime([uxind],tdata.time)[0]
      tdataDat = np.concatenate(tdata._getStepsShots(tstride[0],tstride[1]),axis=0)
      repack = unravelScanSteps(repack, [len(txind) for txind in xind])
      package.append((tdataDat,repack))

    #unpack single event steps and calculate
    td = []
    for iNo,tnind in enumerate(nind):
      ipdat = []
      for pack in package:
	tdataDat,repack = pack
	ipdat.append(tdataDat[tools.smartIdx(repack[iNo])])
      
      td.append(func(ipdat))

    output.append((np.asarray(td),))
  return output


    






      





   


  pass

def get_common_timestamps(allobjects):
  times = None
  for o in allobjects:
    if times==None:
      times = o.time
    else:
      ia,io = filterTimestamps(times,o.time)
      iar,stpsz = ravelScanSteps(ia)
      times = np.hstack([tt for tt in times if len(tt)>0])[iar]
      times = unravelScanSteps(times,stpsz)
  return times

############## HELPER #########################
class Evaluate(object):
  def __init__(self,datainstance):
    self.data = datainstance

  def _checkParent(self):
    if self.data._procObj is not None:
      args = self.data._procObj['args']
      for targ in args:
        if isinstance(targ,data):
          continue
        else:
          print "fund evaluating first",tarrg
          targ.evaluate()

  def __call__(self,stepSlice=slice(None),evtSlice=slice(None),progress=True,force=False):
    """Process all data of data instance which will be saved to a ixp hdf5 file. This function could be candidate for parallelization or background processing
    """

    allchunks = self.data._memIterate(stepSlice,evtSlice)
    timestamps = []
    for tstep,tstepchunk in allchunks:
      tts = [self.data.time[tstep][tchunk] for tchunk in tstepchunk if len(tchunk)>0]
      if len(tts)>0:
	timestamps.append(np.concatenate(tts))
      else:
	timestamps.append([])
    #timestamps = [np.concatenate([self.data.time[tstep][tchunk] for tchunk in tstepchunk if len(tchunk)>0]) for tstep,tstepchunk in allchunks]

    ixp,path = self.data._getIxpHandle()
    if path in ixp.fileHandle and not force:
      deleteit = 'y' == raw_input("Dataset %s exists in ixp file, would you like to delete it? (y/n) "%path)
      if deleteit:
	del ixp.fileHandle[path]
      else:
	return
    elif path in ixp.fileHandle and force:
      del ixp.fileHandle[path]


    
    grp = ixp.fileHandle.require_group(path)
    grp['_dtype'] = 'data'
    ixp.save(timestamps,grp,name='time')
    dh = grp.require_group('data')
    eventShape = self.data._sizeEvt['shape']
    # progress bar
    Nevtot = np.sum([ np.sum([len(tchunk) for tchunk in step]) for stepNo,step in allchunks])
    processedevents = 0
    widgets = ['Evaluating: ', pb.Percentage(), ' ', pb.Bar(),' ', pb.ETA(),'  ']
    pbar = pb.ProgressBar(widgets=widgets, maxval=Nevtot).start()
    for stepNo,step in allchunks:
      totlen = np.sum([len(x) for x in step])
      totshape = tuple([totlen] + list(eventShape))
      startind = 0
      for chunk in step:
	if len(chunk)==0:
	  tdat = []
	else:
          tdat = self.data[stepNo,np.ix_(chunk)][0]
	if totshape[0]==0:
	  print stepNo,totshape
          ixp.save([],grp,name='time/data/#%06d'%stepNo)
	else:
	  if len(tdat)>0:
	    ds = grp.require_dataset('data/#%06d'%stepNo,totshape,dtype=tdat.dtype)
	    ds[startind:startind+len(chunk),...] = tdat
        processedevents += len(chunk)
        pbar.update(processedevents)

	ixp.fileHandle.flush()

	startind += len(chunk)
    pbar.finish()
    
    #raise NotImplementedError('Use the source, luke!')
    #self.data = Data(time=timestamps,ixpAddressData=grp['data'],name=self.data.name)
    self.data._ixpAddressData = grp['data'] 
    self.data._rdStride = self.data._rdFromIxp
    self.data._procObj = None
    self.data._time = timestamps
    lens = [len(td) for td in self.data._time]
    self.data._filter = unravelScanSteps(np.arange(np.sum(lens)),lens)
    self.data._Nsteps = len(lens)
    self.data._lens = lens

  def __getitem__(self,x):
    n = self.data._Nsteps
    #lens = [len(tfilt) for tfilt in self._filter]
    x = tools.iterfy(x)
    if len(x)==1:
      if n==1:
	evtSlice = x[0]
	stepSlice = slice(None)
      else:
	evtSlice = slice(None)
	stepSlice = x[0]
    elif len(x)==2:
      evtSlice = x[1]
      stepSlice = x[0]
    
    return self.__call__(stepSlice=stepSlice,evtSlice=evtSlice)

  ## TODO: function that evaluates random sample of events...

class scanVar(object):
  def __init__(self,fhandle,name,paths):
    self._h5s = fhandle
    self._name = name
    self._paths = paths
    self._checkNcalib()
    self._read()

  def _read(self):
    names = self._h5s[0][self._paths % 0]["name"]
    self.name  = list(names)[0]
    self.names = list(names)
    for i in range(len(self._h5s)):
      h=self._h5s[i]
      v = []
      for c in range(self._numOfScanSteps[i]):
        val = h[self._paths % c]["value"]
        for (n,v) in zip(names,val):
          tools.addToObj(self,address(i,c,n),v)
    for n in names:
      dhelp = memdata(self,n,0)
      tools.addToObj(self,n,dhelp)

  def __repr__(self):
    return "`scanVar` object: %s" % self._name

  # data object (self)manipulation

  def _checkNcalib(self):
    self._numOfScanSteps = []
    for h in self._h5s:
      n=0
      while( tH5.datasetExists(h,self._paths % n) ):
        n+=1
      self._numOfScanSteps.append(n)
    return self._numOfScanSteps

  def __getitem__(self,x):
    return tools.getFromObj(self,self.name)[x]

def address(fileNum,stepNum,what):
  return "_file%d.step%d.%s" % (fileNum,stepNum,what)

def initmemdataraw(data):
  if data==None:
    return None,None
  if type(data) is dict:
    dk = data.keys()
    if not (('data' in dk) and ('time' in dk)):
      print "Raw memdata should be dict with keys data and time of same length"
      return
    dat = data['data']
    tim = data['time']
  elif type(data) is list:
    dat = data[0]
    tim = data[1]
  if tim==None:

    print "NB: memdata instance without timestamps!"
  elif not len(dat)==len(tim):
    print "data and time in raw memdata should have same length"
    return
  return dat,tim

def ravelScanSteps(ip):
  stepsizes = [len(tip) for tip in ip]
  rip = [tip for tip in ip if len(tip)>0]
  if len(rip)==0:
    op = []
  else:
    op = np.hstack(rip)
  return op,stepsizes

def ravelIndexScanSteps(ip,stepsizes,stepNo=None):
  if (not len(ip)==len(stepsizes)):
    if not stepNo==None:
      stepNo = tools.iterfy(stepNo)
      aip = [[]]*len(stepsizes)
      for tstepNo,tip in zip(stepNo,ip): 
	aip[tstepNo] = tip
      ip = aip
    else:
      raise Exception('index list doesn\'t fit stepsizes; Extra keyword argument stepNo required.')
  op = []
  asz= 0
  for tip,tsz in zip(ip,stepsizes):
    tip = np.asarray(tip,dtype=int)
    op.append(tip+asz)
    asz+=tsz
  return np.hstack(op)

def unravelScanSteps(ip,stepsizes):
  op = []
  startind = 0
  stopind = 0
  for stepNO,stepsize in enumerate(stepsizes):
    stopind += stepsize
    op.append(ip[startind:stopind])
    startind += stepsize
  return op

def unravelIndexScanSteps(ip,stepsizes,replacements=None):
  op = []
  startind = 0
  stopind = 0
  for stepsize in stepsizes:
    stopind += stepsize
    tind = np.logical_and(startind<=ip,ip<stopind)
    if replacements==None:
      op.append(ip[tind])
    else:
      op.append(replacements[tind])
    startind += stepsize
  return op

def getStepShotsFromIndsTime(inds,times,stride=None,getIncludedSteps=False):
  """ Takes inds (indices of the ravelled output sorted in lists of steps) as well as the time array of a detector to sort into inds. returns list for each nonempty inds step sontaining array of involved steps and shots in those steps.
  """
  addresses = [[(step,shot) for shot in range(len(ttime))] for step,ttime in enumerate(times)]
  added = []
  for addr in addresses:
    added.extend(addr)
  added = np.asarray(added)
  stepShots = []
  #print inds
  includedsteps = []
  for step,ind in enumerate(inds):
    if len(ind)>0:
      includedsteps.append(step)
      ts = np.vstack(added[np.ix_(ind)])
      if not stride==None:
	if len(stride)>step and np.iterable(stride[step]):
	  tstride = stride[step]
	else:
	  tstride = stride
	if np.max(tstride)<len(ts):
	  ts = ts[tstride]
      tsteps = np.unique(ts[:,0])
      tshots = [ts[(ts[:,0]==n).nonzero()[0],1] for n in tsteps]
      stepShots.append([tsteps,tshots])

  #print stepShots
  if getIncludedSteps:
    return stepShots,np.asarray(includedsteps)
    
  return stepShots

def matchEvents(*args):
  filt = args[0].ones()
  for tobj in args[1:]:
    filt *=tobj.ones()
  ret = [filt*tobj for tobj in args]
  return tuple(ret)

def getClosestEvents(ts0,ts1,N,excludeFurtherThan=None):
  """finds N closes events in ts0 to ts1. Returns those events (after 
  ravelling) plus the unravel information (lengths per step in ts0 events).
  """
  ts0r,ss0 = ravelScanSteps(ts0)
  ts1r,ss1 = ravelScanSteps(ts1)
  ts0r = getTime(ts0r,asTime=True)
  ts1r = getTime(ts1r,asTime=True)
  indout = [np.sort(np.abs(ts0r-tts1r).argsort()[:N]) for tts1r in ts1r]
  return indout,ss0

def getTime(tsi,resolution = 1e-9,asTime=False):
  """Takes structured second/nanosecond array and returns integer-clipped ms array.
  """
  if tsi.dtype.names:
    nam = tsi.dtype.names
    if (not asTime) and (('fiducials' in nam) and ('seconds' in nam)):
      tso = tsi
    elif ('nanoseconds' in nam) and ('seconds' in nam):
      tso = np.int64(tsi['seconds'])*int(1/resolution) + tsi['nanoseconds']/int(1e-9/resolution)
  else:
    tso = tsi
  return tso

def filterTimestamps(ts0,ts1,asTime=True):
  """returns 2 list of indices arrays: 
    1) ts0 subset that is in ts1 
    2) bring ts1 in the order of ts0"""
  #raise NotImplementedError('Use the source, luke!')
  ts0r,ss0 = ravelScanSteps(ts0)
  ts1r,ss1 = ravelScanSteps(ts1)
  ts0r = getTime(ts0r,asTime=asTime)
  ts1r = getTime(ts1r,asTime=asTime)
  
  if (not (ts0r.dtype.names==None or ts1r.dtype.names==None)) and not asTime:
    sel0rb,sel1ri = sortTimestampFiducials(ts0r,ts1r)
  else:
    sel0rb = np.in1d(ts0r,ts1r)
    sel1rb = np.in1d(ts1r,ts0r)
    sel1ri = sel1rb.nonzero()[0][ts1r[sel1rb].argsort()[ts0r[sel0rb].argsort().argsort()]]
  op0 = unravelIndexScanSteps(sel0rb.nonzero()[0],ss0)
  ss1a = [len(top0) for top0 in op0]
  op1 = unravelScanSteps(sel1ri,ss1a)
  return op0,op1

def sortTimestampFiducials(a,b,P=300,maxoff=10,remove_duplicates=False,remove_all_doubles=False):
  """ sorting based on timestamps AND fiducials
  brings those b that are in a into order of a
  Takes the fiducial period (in ts units,<=the real period) and the max tx 
  offset betwen a and b. Returns a boolead array of what in 
  a is also in b and an array of indices that bring b in the 
  order of a (only for those that exist there)"""

  ats,af = (a['seconds'],a['fiducials'])
  bts,bf = (b['seconds'],b['fiducials'])
  # Getting timestamp groups
  # down to integer units of the smallest offset to be considered
  atsmn = min(ats)
  ared = np.floor((ats-atsmn)/maxoff)
  bred = np.floor((bts-atsmn)/maxoff)
  #get number of small groups per search group (a), as well 
  #as the number of search groups.
  Na = np.floor(P/maxoff)-2 # number of small groups per search group in a
  Nintervals = np.ceil(max(ared)/Na) # number of search groups
  av = np.arange(Nintervals)  #this can be used for index comparison when looping through search groups
  ai = np.searchsorted(np.arange(Nintervals)*Na , ared ,side='right')-1
  print Na,Nintervals
  print np.floor((ats-atsmn)/maxoff)
  print av
  print ai
  Sela = []
  Selbri = []
  for nInt in av:
    aInd = ai==nInt
    bInd = np.logical_and(bred >= (nInt*Na-1), bred <= ((nInt+1)*Na+1))
    taf = af[aInd]
    tbf = bf[bInd]
    if remove_duplicates:
      aindDup = removeDuplicates(taf,remove_all_doubles=remove_all_doubles)
      bindDup = removeDuplicates(tbf,remove_all_doubles=remove_all_doubles)
    else:
      aindDup = np.ones(np.shape(taf),dtype=bool)
      bindDup = np.ones(np.shape(tbf),dtype=bool)
      print 'notlookingfordoubles'
    taf = taf[aindDup]
    tbf = tbf[bindDup]
    if not len(taf)==len(np.unique(taf)):
      print "doubles in taf"
    if not len(tbf)==len(np.unique(tbf)):
      print "doubles in tbf"

    sela = np.in1d(taf,tbf)
    selb = np.in1d(tbf,taf)
    selbri = selb.nonzero()[0]\
	[tbf[selb].argsort()[taf[sela].argsort().argsort()]]
    selbri = bInd.nonzero()[0][bindDup][selbri]
    aindDup[aindDup] = sela
    Sela.append(aindDup)
    print len(aindDup)
    Selbri.append(selbri)
  return np.hstack(Sela),np.hstack(Selbri)

def removeDuplicates(a,remove_all_doubles=True):
  asrti = a.argsort()
  asrt  = a[asrti]
  cleansorted = np.concatenate(([True], asrt[1:] != asrt[:-1]))
  if remove_all_doubles:
    cleansorted = np.concatenate((cleansorted, [True]))
    cleansorted = np.logical_and(cleansorted[:-1],cleansorted[1:])
  return cleansorted[asrti.argsort()]

def getScanVec(instance):
  leninst = len(instance)
  scan = instance.scan
  names = np.asarray(scan._get_keys())
  lens = np.asarray([len(scan[tname]) for tname in names])
  isbincenter = np.asarray([tname.split('_')[-1]=='bincenter' for tname in names])
  if sum(lens==leninst)<1:
    print "Attention: no suitable scan vector found for %s !" %instance._name
    #raise Exception
    return None
  elif sum(isbincenter)==1:
    sel_name = names[isbincenter][0]
    return sel_name, scan[sel_name]
  elif sum(isbincenter)>1:
    sel_name = names[isbincenter[0]][0]
    print "Attention: More than one bincenter array found for %s, will use %s." \
	%(instance._name,sel_name)
    return sel_name, scan[sel_name]
  else:
    sel_name = names[lens==leninst][0]
    return sel_name, scan[sel_name]
  
def digitizeN(*args,**kwargs):
  """multi-dimensional digitization, used memdata or timestamps to sort data"""
  target = kwargs.get('target', None)
  if isinstance(target,memdata):
    dat,stsz = ravelScanSteps(target.data)
    tim,stszt = ravelScanSteps(target.time)
  elif target==None:
    atarg = args[0].ones()
    for targ in args[1:]:
      atarg = atarg*targ.ones()
    dat,stsz = ravelScanSteps(atarg.data)
    tim,stsz = ravelScanSteps(atarg.time)
  else: #assuming now it is timestamps
    tim,stsz = ravelScanSteps(target)
    dat = np.ones(np.shape(tim))
  indmat = np.nan*np.ones((len(tim),len(args)))
  vecs = []
  names = []

  for narg,targ in enumerate(args):
    tname,tvec = getScanVec(targ)
    vecs.append(tvec)
    names.append(tname)
    isused,sortThis = filterTimestamps(targ.time,[tim],asTime=True)
    for nstep,sortThisStep in enumerate(sortThis):
      indmat[sortThisStep,narg] = nstep
  gd = ~np.isnan(indmat).any(axis=1)
  indmat = indmat[gd,:]
  #print np.shape(indmat)
  #from matplotlib import pyplot as plt
  #plt.clf()

  totshape = tuple([len(tvec) for tvec in vecs])
  indmat = np.asarray(indmat.T,dtype=int)
  grouping = np.ravel_multi_index(indmat,totshape)
  groupingall = np.ones(len(tim),dtype=int)
  groupingall[gd] = grouping

  timout = [ tim[groupingall==i] for i in range(np.prod(totshape))]
  
  datout = [ dat[groupingall==i] for i in range(np.prod(totshape))]
  scan = tools.dropObject(name='scan')

  meshes = tools.ndmesh(*tuple(vecs))
  for tname,tmesh in zip(names,meshes):
    scan[tname] = tmesh.ravel()

  grid = Grid(zip(names,vecs))


  return memdata(name='%d-dimensional histogram'%(len(totshape)),input=[datout,timout],scan=scan,grid=grid)
      
class Grid(object):
  def __init__(self,definition):
    self._definition = definition
    self._ixpsaved = ['_definition']

  def _get_vec(self,sel=None):
    if sel==None:
      return tuple([tdef[1] for tdef in self._definition])
    elif type(sel)==int:
      return self._definition[sel][1]
    elif type(sel)==str:
      ind = (sel==np.asaray(self.names)).nonzero()
      return self._definition[ind][1]

  def _get_name(self,sel=None):
    if sel==None:
      return tuple([tdef[0] for tdef in self._definition])
    elif type(sel)==int:
      return self._definition[sel][0]
  names = property(_get_name)

  def _get_shape(self):
    vecs = self._get_vec()
    return tuple([len(tv) for tv in vecs])
  shape = property(_get_shape)    
 
  def gridshape(self,vec): 
    if not len(vec)==np.prod(self.shape):
      raise Exception('length of vec does not match grid size!')
    else:
      vec = np.asarray(vec)
      return np.reshape(vec,self.shape)


  def format(self,vec,dims=None,method=np.mean):
    isSelection = np.iterable(vec[0])
    vec = self.gridshape(vec)
    if not dims==None:
      adims = list(self.shape)
      for rdim in dims:
	adims.remove(rdim)
      trsp  = dims+tuple(adims)
      vec = vec.transpose(trsp)
    return self[:]+(vec,)
 
      #if len(adims)>0:
	#stridestr = '[...'+len(adims)*',:'+']'
	#if isSelection:
	  #exec('vec'+stridestr+' = np.concatenate(vec'+stridestr+'.ravel(),axis=0)')
	#else:
	  #vec[len(dims)-1] = np.concatenate(vec[len(dims)-1])










  def __len__(self):
    return len(self._definition)

  def __getitem__(self,x):
    if type(x)==slice:
      inds = range(*x.indices(len(self)))
      return tuple([self._get_vec(sel=ind) for ind in inds])
    else:
      return self._get_vec(sel=x)

############ IXP cache file functions #############
_ixpOpenHandles = []
class Ixp(object):
  def __init__(self,filename):
    self.fileName = filename
    if os.path.isfile(self.fileName):
      print 'Found cached data in %s' %(self.fileName)
    self._fileHandle = None
    self._forceOverwrite = False

  def get_cacheFileHandle(self,reopen=False,mode='a'):
    if self._fileHandle is None:
      if reopen:
	openNames = [th.filename for th in _ixpOpenHandles]
	if self.fileName in openNames:
	  exHandle = _ixpOpenHandles[openNames.index(self.fileName)]
	  exHandle.close()
	  _ixpOpenHandles.remove(exHandle)
      cfh = h5py.File(self.fileName,mode)
      self._fileHandle = cfh
      _ixpOpenHandles.append(cfh)
    else:
      cfh = self._fileHandle
      if reopen:
	if cfh in _ixpOpenHandles:
	  sfh.close()
	  _ixpOpenHandles.remove(exHandle)
          self._fileHandle = h5py.File(self.fileName,mode)
      cfh = self._fileHandle
    return cfh

  def set_cacheFileHandle(self,cfh):
    self._fileHandle = cfh

  fileHandle = property(get_cacheFileHandle,set_cacheFileHandle)
 
  def delete(self):
    ristr ='Do you intend to permanently delete file \n  %s\n  (y/n) ' %(self.fileName)
    if 'y'==raw_input(ristr):
      os.remove(self.fileName)

  #def save(self, obj, filename='default', force=False):
    #if not filename=='default':
      #self.fileName = filename
    
    #cfh = self.fileHandle
    ## datasets in "root"-dataset instance
    #fGroup = cfh.require_group('base_instance')
    #for sfield in obj._ixpsaved:
      #if sfield in fGroup.keys() and not force:
        #rawstr = 'Overwrite %s in base dataset ? (y/n/a) ' %(sfield)
        #ret = raw_input(rawstr)
        #if ret=='a': del fGroup[sfield]; force = True
        #if ret=='y': del fGroup[sfield]
        #if ret=='n': continue
      #elif sfield in fGroup.keys() and force:
        #print "about to delete %s" %(sfield)
        #del fGroup[sfield]
      #self.mkDset(fGroup,sfield,self.config.base.__dict__[sfield])

    ## detectors and data datasets
    #for field in obj._ixpsaved:
      #fGroup = cfh.require_group(field)
      #for sfield in self.config.base.__dict__[field]._savelist:
        #if sfield in fGroup.keys() and not force:
          #rawstr = 'Overwrite %s in %s ? (y/n/a) ' %(sfield, field)
          #ret = raw_input(rawstr)
          #if ret=='a': del fGroup[sfield]; force = True
          #if ret=='y': del fGroup[sfield]
          #if ret=='n': continue
        #elif sfield in fGroup.keys() and force:
          #print "about to delete %s" %(sfield)
          #del fGroup[sfield]
        #self.mkDset(fGroup,sfield,self.config.base.__dict__[field].__dict__[sfield])
    #cfh.close()
  
  def save(self, obj, parenth5handle, name=None, force=None):
    if hasattr(obj,'_isIxp') and obj._isIxp():
      return
    if force is None:
      force = self._forceOverwrite
    else:
      self._forceOverwrite = force
    pH = parenth5handle
    isgr = hasattr(obj,'_ixpsaved')
    if name is None:
      name = obj._name
    if isgr:
      fGroup = pH.require_group(name)
      for sfield in obj._ixpsaved:
	if not hasattr(obj,sfield[0]):
	  print "Error: field %s in group %s is not existing!" %(sfield,name)
	  continue
        if not sfield[1]=='no':
          self.save(obj.__dict__[sfield[0]],fGroup,name=sfield[0])
    else:
      if name in pH.keys() and not force:
        rawstr = 'Overwrite %s in %s ? (y/n/a) ' %(name,pH.name)
        ret = raw_input(rawstr)
        if ret=='a': del pH[name]; force = True; self._forceOverwrite = True
        if ret=='y': del pH[name]
        if ret=='n': pass
      elif name in pH.keys() and force:
        print "about to delete %s" %(name)
        del pH[name]
      self.mkDset(pH,name,obj)


  def load(self,parenth5handle=None):
    #savobjs = self.cacheFileHandle.keys()
    pH = parenth5handle
    if pH is None:
      pH = self.fileHandle
    isgr = isinstance(pH,h5py.Group) and (not '_dtype' in pH.keys())
    name = os.path.split(pH.name)[-1]
    if isgr:
      op = tools.dropObject(name=name)
      for field in pH.keys():
	ob = self.load(pH[field])
	op.__setitem__(field,ob,setParent=True)
    else:
      op = self.rdDset(pH)
    return op


  def mkH5list(self,listh):
    H5list = []
    keylist = listh.keys()
    keylist = np.sort(keylist)
    for h in keylist:
      H5list.append(listh[h])
    return H5list

  def mkDset(self,rootgroup,name,data):
    #if name=='scanMot':
      #de=bug
    if isinstance(data,memdata):
      if hasattr(data,'_data') and hasattr(data,'data') and hasattr(data,'time'):
	newgroup = rootgroup.require_group(name)
	try:
	  newgroup.create_dataset('_dtype',data='memdata')
	  self.mkDset(newgroup,'data',data.data)
	  self.mkDset(newgroup,'time',data.time)
	  self.save(data.scan,newgroup,name='scan')
	  if not data.grid==None:
	    gridgroup = newgroup.require_group('grid')
	    self.mkDset(gridgroup,'definition',data.grid._definition)
	  
	except Exception,e:
	  print e
	  print "could not save memdata instance %s" %name
	  del newgroup

    if type(data) is list:
      newgroup = rootgroup.require_group(name)
      newgroup.create_dataset('_dtype',data='list')
      el=0
      for d in data:
        tname = '#'+'%06d' %(el)
        self.mkDset(newgroup,tname,d)
        el += 1
    if type(data) is tuple:
      newgroup = rootgroup.require_group(name)
      newgroup.create_dataset('_dtype',data='tuple')
      el = 0
      for d in data:
        tname = '#%06d' %el
        self.mkDset(newgroup,tname,d)
        el += 1
    if type(data) is dict:
      newgroup = rootgroup.require_group(name)
      newgroup.create_dataset('_dtype',data='dict')
      for key in data.keys():
        tname = '%s' %key
        self.mkDset(newgroup,tname,data[key])
    if type(data) in [np.ndarray,int,float,str,np.unicode_,
                      np.string_,np.bool,np.int,
                      np.int8,np.int16,np.int32,np.int64,
                      np.uint8,np.uint16,np.uint32,np.uint64,
                      np.float,np.float32,np.float64,
                      np.complex,np.complex64,np.complex128,
                      np.core.records.recarray	
                      ]:
      #if type(data) is str: print data; print rootgroup.name
      #newgroup = rootgroup.require_group(name)
      #newgroup.create_dataset('_dtype',data='')

      if not ((data==[]) or (type(data)==np.ndarray and data.size==False)):
        rootgroup.create_dataset(name,data=data)
      else:
        rootgroup.create_dataset(name,data='empty')

  def rdDset(self,rootgroup):
    try:
      allkeys = rootgroup.keys()
    except:
      allkeys = []
    if '_dtype' in allkeys:
      tdtype = rootgroup['_dtype'].value
      if tdtype == 'memdata':
	name = os.path.split(rootgroup.name)[-1]
	try:
	  data = self.rdDset(rootgroup['data'])
	  time = self.rdDset(rootgroup['time'])
	  try:
	    scan = self.load(parenth5handle=rootgroup['scan'])
	  except Exception,e:
	    print e
	    scan = None
	  if 'grid' in rootgroup.keys():
	    try:
	      grid = Grid(self.rdDset(rootgroup['grid']['definition']))
	    except Exception,e:
	      print e
	      grid = None
	  else:
	    grid = None
	  return memdata(input = [data,time], name=name, scan=scan, grid=grid)
	except Exception,e:
	  print e
	  print "reading error with memdata instance %s" %name
	  return name+"_empty_corrupt"
      if tdtype == 'data':
	name = os.path.split(rootgroup.name)[-1]
	try:
	  ixpAddressData = rootgroup['data']
	  time = self.rdDset(rootgroup['time'])
	  return Data(time=time,ixpAddressData=ixpAddressData,name=name)
	except Exception,e:
	  print e
	  print "reading error with memdata instance %s" %name
	  return name+"_empty_corrupt"
      if tdtype == 'list':
        listkeys = [key for key in allkeys if '#' in key]
        listkeys.sort()
        data = []
        for key in listkeys:
	  #data.append(self.rdDset(rootgroup[key]))
          keydat = self.rdDset(rootgroup[key])
          data.append(keydat)

      if tdtype == 'tuple':
        listkeys = [key for key in allkeys if '#' in key]
        listkeys.sort()
        data = []
        for key in listkeys:
	  #data.append(self.rdDset(rootgroup[key]))
          keydat = self.rdDset(rootgroup[key])
          data.append(keydat)
	data = tuple(data)
      
      if tdtype == 'dict':
        dictkeys = [key for key in allkeys if '_dtype' not in key]
        data = dict()
        for key in dictkeys:
          data[key] = self.rdDset(rootgroup[key])

    else:
      data = rootgroup.value
      if data=='empty':
        data = []

    return data

  def appendData(self,datagroup,data,step=None):
    assert len(data)==2
    tim = data[1]
    dat = data[0]
    shp = dat.shape
    try:
      datagroup['_dtype'] = 'data'
    except:
      del datagroup['_dtype']
      datagroup['_dtype'] = 'data'
    dh = datagroup.require_group('data')
    dt = datagroup.require_group('time')
    try:
      datagroup['time/_dtype'] = 'list'
    except:
      del datagroup['time/_dtype']
      datagroup['time/_dtype'] = 'list'
    dhk = [int(tmp[1:]) for tmp in dh.keys() if tmp[0]=='#']
    dtk = [int(tmp[1:]) for tmp in dt.keys() if tmp[0]=='#']
    assert dhk == dtk, "Inconsistency in data instance ixp file!"
    if step is None:
      if len(dhk)==0:
	step = 0
      else:
        step = max(dhk)+1

    if step in dtk:
      dsd = datagroup['data/#%06d'%step]
      dst = datagroup['time/#%06d'%step]
      assert dsd.shape[0]==dst.shape[0]
      assert dsd.shape[1:]==shp[1:]
      dsd.resize((dsd.shape[0]+shp[0],)+shp[1:])
      dst.resize((dst.shape[0]+shp[0],))
      dsd[-shp[0]:] = dat
      dst[-shp[0]:] = tim

    else:
      dsd = datagroup.create_dataset('data/#%06d'%step,data=dat,maxshape=(None,)+shp[1:])
      dst = datagroup.create_dataset('time/#%06d'%step,data=tim,maxshape=(None,))

class mem(object):
  def __init__(self):
    self.f=open("/proc/meminfo","r")
    tot = self.f.readline()
    self.tot = self.stringToBits(tot)
    self.updatefree()

  def updatefree(self):
    self.f.seek(0); # rewind
    tmp = self.f.readline()
    memfree = self.f.readline()
    return self.stringToBits(memfree)

  def stringToBits(self,s):
    return 1024*8*int(s.split()[1])

  free = property(updatefree)

############ CONFIG FILE #############
def _rdConfigurationRaw(fina="ixppy_config"):
  """Reads configuration file that has the description of the different detectors, the default data paths etc."""
  if not os.path.isfile(fina):
    path = os.path.abspath(__file__)
    fina = os.path.dirname(path) + '/' + fina    
  filecontent = dict()
  fhandle = open(fina)
  lines = fhandle.readlines()
  fhandle.close()
  foundlabel = False
  lines.append("#* end"); # needed to read the last block
  for i in range(len(lines)):
    line = lines[i].strip()

    if foundlabel and (len(line)!=0):
      if line[0] is not '#':
        if line.split():
          tdat.append(line.split())
      else:
        # case when dataset should be finished
        filecontent[tlabel] = tdat

    #label line, keyword to characterize dataset starting next line
    if line[:2]=='#*':
      line = line[2:]
      tlabel = line.split()[0]
      foundlabel = True
      tdat = []
  return filecontent

def rdConfiguration(fina="ixppy_config",beamline=None):
  dat = _rdConfigurationRaw(fina)
  home = expanduser("~")
  if os.path.exists(home+'/.ixppyrc'):
    datuser = _rdConfigurationRaw(home+'/.ixppyrc')
    dat = tools.dictMerge(dat,datuser)
  cnf = dict()
  # structurize file output to enable later restructuring of configuration file
  # dataset aliases

  # read present beamline
  if beamline==None:
    cnf["beamline"] = dat["beamline"][0][0]
  else:
    cnf["beamline"] = beamline

  pointDetbeamline = _interpreteDetectorConfig(dat,'pointDet',cnf['beamline'])
  pointDetcommon   = _interpreteDetectorConfig(dat,'pointDet','common')
  pointDetspecial   = _interpreteDetectorConfig(dat,'pointDet','special')
  cnf["pointDet"] = dict(pointDetcommon.items() + pointDetbeamline.items() + pointDetspecial.items())
  areaDetbeamline = _interpreteDetectorConfig(dat,'areaDet',cnf['beamline'])
  areaDetcommon   = _interpreteDetectorConfig(dat,'areaDet','common')
  cnf["areaDet"] = dict(areaDetcommon.items() + areaDetbeamline.items())

  #cnf["areaDet"] = _interpreteDetectorConfig(dat,'areaDet',cnf['beamline'])
  cnf["dataPath"] = _interpreteCnfPath(dat,cnf['beamline'],what="data_path")
  cnf["cachePath"] = _interpreteCnfPath(dat,cnf['beamline'],what="cache_path")
  cnf["scan_step"] = dat['scan_step']
  cnf["cache_directory"] = dat['cache_directory'][0]
  cnf["epics_dset"] = dat['epics_dset'][0]
  cnf["epics_cc_dset"] = dat['epics_dset'][0][1]

  return cnf

############ FILE/PATH TOOLS  #############
def _interpreteCnfPath(confraw,beamline,what="data_path"):
  cnf = dict()
  if (what == "cache_path"):
    for line in confraw["cache_directory"]:
      cnf[line[1]] = line[2]
  else:
    for line in confraw["default_datapath_hosts"]:
      if line[0]==beamline:
        cnf[line[1]] = line[2]
    for line in confraw["default_datapath"]:
      if line[0]==beamline:
        cnf['default'] = line[1]
        break
  return cnf

def _getInputDirectory(cnfFile):
    knownhosts = cnfFile['dataPath'].keys()
    if gethostname() in knownhosts:
      inputDirectory = cnfFile['dataPath'][gethostname()]
    else:
      inputDirectory = os.path.join(cnfFile['dataPath']['default'],cnfFile['beamline'])
    return inputDirectory

def _interpreteDetectorConfig(confraw,detkind,beamline):
  """interprets detector configuration"""
  dat = confraw
  detname = '%s_%s' %(detkind,beamline)
  try:
    cnf = dict()
    aliases = [tdat[0] for tdat in dat[detname]]
    #print aliases
    n=0
    noofbak = dict()
    for alias in aliases: 
      if alias in cnf.keys():
        if alias not in noofbak.keys(): noofbak[alias] = -1
        noofbak[alias] += 1 
        alias = alias + '_bak' + str(noofbak[alias])
      cnf[alias] = dict()
      cnf[alias]['data'] = dat[detname][n][1]
      cnf[alias]['timestamp'] = dat[detname][n][2]
      cnf[alias]['conf'] = dat[detname][n][3]
      n+=1
    return cnf
  except:
    #print "No %s detector datasets found in configuration file!" %(detname)
    return dict()


def _rdExperimentList(fina="xpp_experiment_list.txt"):
  """Reads experiment list with pi names, experiment numbers and experiment titles"""
  if ~os.path.isfile(fina):
    path = os.path.abspath(__file__)
    fina = os.path.dirname(path) + '/' + fina    
  d = np.loadtxt(fina,dtype=str,delimiter='\t')
  l = dict(no=d[:,0],id=d[:,1],startdate=d[:,2],enddate=d[:,3],pi=d[:,4],title=d[:,5])
  return l

def getExperimentnumber(searchstring,newest=True,searchtitle=False):
  l = _rdExperimentList()
  lines = []
  n=0
  for name in l["pi"]:
    if searchstring.lower() in name.lower():
      lines.append(n)
    n+=1
  if len(lines)>1:
    sd = []
    for tsd in l["startdate"][np.ix_(lines)]:
      sd.append(dateutil.parser.parse(tsd))
    sind = np.argsort(sd)
    pi = l["pi"][lines[0]]
    id = l["id"][np.ix_(lines)][sind]
    no = l["no"][np.ix_(lines)][sind]
    ti = l["title"][np.ix_(lines)][sind]
    st = l["startdate"][np.ix_(lines)][sind]
    str =  "Found %d experiments for PI (%s):" %(len(lines),pi)
    for n in range(len(lines)):
      str += "\n  %s\t%s" %(no[n],ti[n])
    print str
  else:
    id = l["id"][np.ix_(lines)]
    no = l["no"][np.ix_(lines)]
  return list(no),list(id)

def getExperimentList(printit=True,sortDates=True):

  l = _rdExperimentList()
  if sortDates:
    ind = np.argsort(l['startdate'])

  if printit:
    for n in ind:
      print "%s    %s    %s"  %(l['startdate'][n],l['no'][n],l['pi'][n])
  return l

############## UNDER CONSTRUCTION ###############

class interp(object):
  ''' Under construction, not working'''
  def __init__(self,memdinst):
    self.md = memdinst
    self._apply_method = None
  def closest(self):
    self._apply_method = 'closest'
  def linear(self):
    self._apply_method = 'linear'
  def closestBefore(self):
    self._apply_method = 'closestBefore'
  def spline(self):
    self._apply_method = 'spline'
  def linearAll(self):
    self._apply_method = 'linearAll'

  def _linear(self,timeOther):
    self._apply_method = None
  def _linearAll(self,timeOther):
    self._apply_method = None
  def _closest(self,timeOther):
    self._apply_method = None
  def _closestBefore(self,timeOther):
    self._apply_method = None
#def dataInterpolate(source,steplengths,sourceInd,stride):
  #for 
  

############## ixppy 1.0 remains ###############

###### EASY filtering tools ##########

class filtDsetObj(object):
  def __init__(self,dset):
    self._dset = dset

  def __getattr__(self,inp):
    self.__dir__()
    return self.__getattribute(inp)
  
  def __dir__(self):
    for det in self._dset.detectors:
      td = self._dset.__dict__[det]
      if hasattr(td,'_masked_arrays'):
	if 'data' in td._masked_arrays.keys():
	  if not td._masked_arrays['data'] == []:
            self.__dict__[det] = filtDetObj(self._dset,det)
    return self.__dict__.keys()

class filtDetObj(object):
  def __init__(self,dset,det):
    self._dset = dset
    self._det = det
    self._channels = self._dset.__dict__[det]._masked_arrays['data']
  
  def __getattr__(self,inp):
    self.__dir__()
    return self.__getattribute(inp)

  def __dir__(self):
    for channel in self._channels:
      channel = channel.strip('_')
      self.__dict__[channel] = filtChannelObj(self._dset,self._det,channel)
    return self.__dict__.keys()

class filtChannelObj(object):
  def __init__(self,dset,det,channel):
    self._dset = dset
    self._det = det
    self._channel = channel
    def channelfun():
      self._dset._initfilter()
      return partial(parameterFilt, self._dset.__dict__[self._det].__dict__['__'+self._channel](),dataset=self._dset)()


    self.__dict__['filt'] = channelfun
    self.__dict__['corr'] = filtchannelcorrobj(self._dset,self._det,self._channel)

class filtchannelcorrobj(object):
  def __init__(self,dset,det,channel):
    self._dset = dset
    self._det = det
    self._channel = channel
  def __getattr__(self,inp):
    self.__dir__()
    return self.__getattribute(inp)

  def __dir__(self):
    self.initalldetchannels()
    return self.__dict__.keys()
  def initalldetchannels(self):
    for det in self._dset.detectors:
      td = self._dset.__dict__[det]
      if hasattr(td,'_masked_arrays'):
	if 'data' in td._masked_arrays.keys():
	  if not td._masked_arrays['data'] == []:
	    self.__dict__[det] = channelcorrobj(self._dset,self._det,self._channel,det)
class channelcorrobj(object):
  def __init__(self,dset,det,channel,corrdet):
    self._dset = dset
    self._det = det
    self._channel = channel
    self._corrdet = corrdet
  def __getattr__(self,inp):
    self.__dir__()
    return self.__getattribute(inp)
  def __dir__(self):
    self._init()
    return self.__dict__.keys()

    

  def _init(self):
    corrdet = self._dset.__dict__[self._corrdet]
    corrchannels = corrdet._masked_arrays['data']
    for cchannel in corrchannels:
      cchannel = cchannel.strip('_')
      def channelfun():
	self._dset._initfilter()
	return partial(parameterCorrFilt,
					  self._dset.__dict__[self._det].__dict__['__'+self._channel](),
					  self._dset.__dict__[self._corrdet].__dict__['__'+cchannel](),
					  dataset=self._dset
					  )()

      self.__dict__[cchannel] = channelfun    

################## data manipulation #################3
def getProfileLimitsNew(Areadet,step=0,shots=range(10),direction=False,lims=None, dname='profile'):
  if lims==None:
    getLims = True
  else:
    getLims = False
  if isinstance(Areadet,data):
    dat = Areadet
    det = None
  else:
    dat = Areadet.data
    det = Areadet
  
  I = np.mean(dat[step,np.ix_(shots)][0],axis=0)
  if direction=='both':
    tools.nfigure('Select limits')
    pl.clf()
    tools.imagesc(I)
    print 'Select region of interest which spans both, horizontal and vertical profile'
    limstot = np.round(tools.getRectangleCoordinates())
    I = tools.subset(I,limstot)

  raise NotImplementedError('Use the source, luke!')
  if direction == 'horizontal' or direction == 'both':
    if getLims:
      tools.nfigure('Select limits')
      pl.clf()
      tools.imagesc(I)
      print 'Select horizontal region of interest'
      lims = np.round(tools.getSpanCoordinates('vertical'))
    if direction == 'both':
      limsver = lims
    else:
      limsdict = dict(projection = direction+' range',limits=lims)
  
  if direction == 'vertical' or direction == 'both':
    if getLims:
      tools.nfigure('Select limits')
      pl.clf()
      tools.imagesc(I)
      print 'Select horizontal region of interest'
      lims = np.round(tools.getSpanCoordinates('horizontal'))
    if direction == 'both':
      limshor = lims
    else:
      limsdict = dict(projection = direction+' range',limits=lims)
  
  if direction == 'both':
    limsdict = dict(projection='box',limits=dict(total=limstot,vertical=limsver,horizontal=limshor))

  tfun = ixppy.wrapFunc(extractProfilesFromData)
  profile = tfun(dat,limsdict)
  if det==None:
    return limsdict,profile
  else:
    if not hasattr(det,dname+'Limits'):
      det._add(dname+'Limits',[])
    det.__dict__[dname+'Limits'].append(limsdict)
    det._add(dname,profile)


def getProfileLimits(Areadet,step=0,shots=range(10),transpose=False,lims=None, dname='profile'):
  if isinstance(Areadet,data):
    dat = Areadet
    det = None
  else:
    dat = Areadet.data
    det = Areadet


  direction = 'vertical'
  if transpose: direction = 'horizontal'
  if lims==None:
    I = dat[step,np.ix_(shots)][0]
    tools.nfigure('Select limits')
    pl.imshow(np.mean(I,axis=0),interpolation='nearest')
    print 'Select region of interest'
    lims = np.round(tools.getSpanCoordinates(direction))
  limsdict = dict(projection = direction+' range',limits=lims)
  tfun = ixppy.wrapFunc(extractProfilesFromData)
  profile = tfun(dat,limsdict)
  if det==None:
    return limsdict,profile
  else:
    if not hasattr(det,dname+'Limits'):
      det._add(dname+'Limits',[])
    det.__dict__[dname+'Limits'].append(limsdict)
    if not direction=='both':
      det[dname] = profile


def extractProfilesFromData(data,profileLimits,cameraoffset=0):
  if profileLimits["projection"]=="vertical range":
    if len(profileLimits['limits'])==0:
      profiles = np.mean(data[:,:,:],axis=1)-float(cameraoffset)
    else:
      profiles = np.mean(data[:,profileLimits['limits'][0]:profileLimits['limits'][1],:],axis=1)-float(cameraoffset)
  elif profileLimits["projection"]=="horizontal range":
    if len(profileLimits['limits'])==0:
      profiles = np.mean(data[:,:,:],axis=2)-float(cameraoffset)
    else:
      profiles = np.mean(data[:,:,profileLimits['limits'][0]:profileLimits['limits'][1]],axis=2)-float(cameraoffset)
  elif profileLimits["projection"]=="box":
    if type(profileLimits['limits'])==dict:
      limstot = profileLimits['limits']['total']
      limshor = profileLimits['limits']['horizontal']
      limsver = profileLimits['limits']['vertical']
      profilehor = np.mean(data[:, 
	min(limsver)+min(limstot[:,1]):max(limsver)+min(limstot[:,1])  ,
	min(limstot[:,0]):max(limstot[:,0])],axis=1)
      profilever = np.mean(data[:  ,
	min(limstot[:,1]):max(limstot[:,1]), 
	min(limshor)+min(limstot[:,0]):max(limshor)+min(limstot[:,0])],axis=2)
      profiles = [profilehor,profilever]

  return profiles

def extractProfilesPeakPar(profiles):
  if type(profiles)==list:
    res = []
    for tprofiles in profiles:
      stack=[]
      for tprofile in tprofiles:
	x = arange(len(tprofile))
	stack.append(tools.peakAna(x,trpofile))
      res.append(np.asarray(stack))
    res = tuple(np.hstack(res).T)
    return res

def extractProfilesPeakParFromData(data,profileLimits):
  return extractProfilesPeakPar(extractProfilesFromData(data,profileLimits))
def getProfilePeakParameters(det):
  data = det.data
  profileLimits = det.profileLimits
  o = ixppy.applyFunction(extractProfilesPeakParFromData,[data,profileLimits],dict(),outputtypes=['memdata']*6,forceCalculation=True)

  det['xpos'] = o[0]
  det['xfwhm'] = o[1]
  det['xpeak'] = o[2]
  det['ypos'] = o[3]
  det['yfwhm'] = o[4]
  det['ypeak'] = o[5]


def rdHdf5dataFromDataSet(h5datasethandle,shots='all'):
  """Reads data from selected shots from dataset handle."""
  dh = h5datasethandle
  #Assume all events/shots are to be read when 'all'
  if type(shots) is str and shots is 'all':
    shots = range(dh.len())
  if type(shots)==np.ndarray:
    shots = list(shots)
  if type(shots) is not list:
    shots = [shots]
  #Make boolean selection array if necessary
  if ~(type(shots)==bool):
    bshots = np.zeros(dh.len(),dtype=bool)
    for i in shots:
      bshots[i] = True
    shots = bshots
  #read data
  dat = dh[shots]
  #if type(dat[1])==tuple:
    #tdat = dict()
    #for tname in dat[1].dtype.names:
      #tdat[tname] = dat[tname]
    #dat=tdat
  return dat


def rdHDFsubset(fina,dataset,idx):

  f = h5py.h5f.open(fina,flags=0)
  dset = h5py.h5d.open(f,dataset)
  dataspace = dset.get_space()
  slabsize = len(idx)
  memspaceID = h5py.h5s.create_simple(tuple([slabsize]),tuple([slabsize]));


  dataspace.select_none()
  for i in idx:
    dataspace.select_hyperslab((i,),1)
  
  datamatrix = dset.read(memspaceID,dataspace)

  dset.close()
  memspaceID.close()
  dataspace.close()
  f.close()

class dropData:
  def add_saved_datafield(self,name,data):
    """adds to a list of datasets to be saved. For faster data recovery and for saving custom analysis progresses (e.g. reduced data from pixel detectors)."""
    if not "_savelist" in self.__dict__:
      self._savelist = []
    self._savelist.append(name)
    self.__dict__[name] = data
  pass

def _timeInMs(time_struc):
  """Makes millisecond array from the read second and nanosecond arrays"""
  ms = np.uint64(time_struc['seconds'])
  ms = ms*1000 + np.round_(time_struc['nanoseconds'],-6)/1000000
  return ms



########## POINT DETECTOR CALCULATIONS ################
def calc_weightedRatio(scanvec,detector,monitor,bins=None,isBinCenter=True):
  if len(scanvec[0])==len(detector[0]):
    #binning
    if not bins:
      bins = np.linspace(np.median(scanvec[0]),
                         np.median(scanvec[-1]),len(scanvec))
    else:
      if len(bins)==1:
        if type(bins) is int:
          bins = np.linspace(np.median(scanvec[0]),
                             np.median(scanvec[-1]),bins)
        elif type(bins) is float:
          bins = np.linspace(np.median(scanvec[0]),
                 np.median(scanvec[-1]),round(len(scanvec)*bins))

    if isBinCenter:
      dbins = np.diff(bins)
      bins = bins + np.hstack([bins[0] - dbins[0]/2., 
                            bins[:-1] + dbins/2., 
                            bins[-1] + dbins[-1]/2.])

    allm = np.hstack(monitor)
    alld = np.hstack(detector)
    alls = np.hstack(scanvec)
    dbin = mbin = wR = []
    inds = digitize(alls,bins[binNo])
    for binNo in range(len(bins)-1):
      tinds = (inds==binNo)
      dbin.append(alld[tinds])
      mbin.append(allm[tinds])
      wR = np.sum(alld[tinds]) / np.sum(allm[tinds])
    return wR,dbin,mbin
  else:
    wR = []
    for td,tm in zip(detector,monitor):
      wR.append(np.sum(td)/np.sum(tm))
  return wR

def digitizeData(data,bins,prcentile=.9):
  """Binning data"""
  dims = len(bins)
  issinglestep = False
  if len(data) is not dims:
    if dims==1:
      data = [data]
    else:

      print "binData: number of data structures does not fit number of bin arrays"
      return
  if iterdepth(data[0])==1:
    ndata = [[td] for td in data]
    data = ndata
    issinglestep = True
    print "in single Step!"

  abins = []
  aN = []
  for dim in range(dims):
    tbins = bins[dim]
    if type(tbins)==int or type(tbins)==float:
      # automatic binning on basis of all cc data
      adata = np.hstack(data[dim])

      if tbins<0:
        # automatic limits
        lims = matplotlib.mlab.prctile(adata,p=((1-prcentile)*100,prcentile*100))
        if type(tbins)==int:
          # number of bins, smart finding of borders, prcentile as kwarg
          cbins = np.linspace(lims[0],lims[1],np.abs(tbins))
        elif type(tbins)==float:
          # bin size, automatic finding of borders, prcentile as kwarg 
          cbins = np.arange(lims[0],lims[1]+np.spacing(1),np.abs(tbins))

      elif tbins==0:
        # graphical input of bin edges
        N,edg = tools.histogramSmart(adata)
        figName = 'Select filter limits'
        tools.nfigure(figName)
        pl.clf()
        pl.step(edg[:-1]+np.diff(edg),N,'k')
        lims = tools.getSpanCoordinates()
        nbins = int(raw_input('Please enter number of bins: '))
        cbins = np.linspace(lims[0],lims[1],np.abs(nbins)+1)

      else:
        # tbins>0
        if type(tbins) is int:
        # number of bins given
          cbins = np.linspace(min(adata),max(adata),tbins+1)

        elif type(tbins) is float:
        # bin size given, full range
          cbins = np.arange(min(adata),max(adata),tbins+1)
        pass

    else:
      # bin edges given
        cbins = tbins

    abins.append(cbins)
    tdata = data[dim]
    tN = []
    for sdata in tdata:
      tN.append(np.digitize(sdata,cbins))
    aN.append(tN)

  # initialize index matrix
  stepbinmatshape = [len(bb)+1 for bb in abins]
  binmatshape = [len(aN[0])]
  binmatshape.extend(stepbinmatshape)
  binmat = np.empty(binmatshape,dtype=np.object_)
  binmat.fill([])
  binmat = np.frompyfunc(list,1,1)(binmat)

  for stepNo in range(len(aN[0])):
    for evNo in range(len(aN[0][stepNo])):
      tind = [stepNo]
      tind.extend([dd[stepNo][evNo] for dd in aN])
      #print tind
      #raw_input()
      binmat[tuple(tind)].append(evNo)
  if issinglestep:
    binmat = binmat[0]

  return binmat,abins

def binStepData(data,binmat):
  bindat = []
  for sdata,sbinmat in zip(data,binmat):
    bindat.append(binData(sdata,sbinmat))
  return np.asarray(bindat)

def binData(data,binmat):
  shp = np.shape(binmat)
  binmat= np.ravel(binmat)
  bindat = []
  for bin in binmat:
    if bin:
      bindat.append(data[np.ix_(bin)])
    else:
      bindat.append([])
  bindat = np.reshape(bindat,shp)
  return bindat

def statStepData(data):
  dmean = []
  dmedian  = []
  dmad  = []
  dstd  = []
  dsum  = []
  dN  = []
  for sdata in data:
    sdmean,sdstd,sdmedian,sdmad,sdsum,sdN = statData(sdata)
    dmean.append(sdmean)
    dstd.append(sdstd)
    dmedian.append(sdmedian)
    dmad.append(sdmad)
    dsum.append(sdsum)
    dN.append(sdN)
    
  return np.asarray(dmean),np.asarray(dstd),np.asarray(dmedian),np.asarray(dmad),np.asarray(dsum),np.asarray(dN)

def statData(data):
  shp = np.shape(data)
  data= np.ravel(data)
  dmean = []
  dmedian  = []
  dmad  = []
  dstd  = []
  dsum  = []
  dN  = []
  for bin in data:
    #stupid numpy workaround...
    try:
      if bin:
        tbe = False
      else:
        tbe=True
    except ValueError:
      if bin.shape==0:
        tbe = True
      else:
        tbe = False

    if not tbe:
      dmean.append(np.mean(bin))
      dmedian.append(np.median(bin))
      dstd.append(np.std(bin))
      dmad.append(tools.mad(bin))
      dsum.append(np.sum(bin))
      dN.append(len(bin))
    else:
      dmean.append(np.nan)
      dmedian.append(np.nan)
      dmad.append(np.nan)
      dstd.append(np.nan)
      dsum.append(np.nan)
      dN.append(np.nan)
  dmean = np.reshape(dmean,shp)
  dmedian = np.reshape(dmedian,shp)
  dmad = np.reshape(dmad,shp)
  dsum = np.reshape(dsum,shp)
  dN = np.reshape(dN,shp)
  dstd = np.reshape(dstd,shp)
  return dmean,dstd,dmedian,dmad,dsum,dN

###### Filtering events ############

def filter(dat,lims=None,graphicalInput=True,figName=None,perc=False):

  if not figName:
    figName = 'Select filter limits'

  if graphicalInput and lims==None:
    tools.nfigure(figName)
    pl.clf()
    N,edg = tools.histogramSmart(dat)
    pl.step(edg[:-1]+np.diff(edg),N,'k')
    lims = tools.getSpanCoordinates()
  elif perc:
    lims = np.percentile(dat,lims)
    
  if type(lims) is bool:
    filt = (dat==lims)
  else:
    filt = (dat>lims[0])&(dat<lims[1])

  return lims,filt

def digitize(dat,bins=None,graphicalInput=True,figName=None):

  if not figName:
    figName = 'Select filter limits'

  if graphicalInput and bins==None:
    tools.nfigure(figName)
    pl.clf()
    N,edg = tools.histogramSmart(dat)
    pl.step(edg[:-1]+np.diff(edg),N,'k')
    lims = tools.getSpanCoordinates()
    ip = 0
    bins = None
    hs = []
    while not ip=='q':
      ip = raw_input('Enter number of bins (q to finish) ')
      if ip=='q': continue
      bins = np.linspace(np.min(lims),np.max(lims),int(ip)+1)
      for th in hs: 
	th.remove()
	hs = []
      hs = [pl.axvline(te) for te in bins]
    print "Selected bins: %g, %g,  "%(np.min(lims),np.max(lims))
  
  return np.digitize(dat,bins),bins

def digitize_new(dat,bins=None,graphicalInput=True,figName=None):

  
  if isinstance(dat,memdata):
    #####
    dat,stsz = ravelScanSteps(dat.data)
    tim,stsz = ravelScanSteps(dat.time)
    #dimension = len(args)
    #res_bins_dat = [dat]
    #res_bins_tim = [tim]
    #res_binvec = []
    #for tbins in args:
      #raise NotImplementeError
      #if isinstance(tbins,memdata):
	#print "Found memdata"
	#res_bins_dat_tmp = []
	#res_bins_tim_tmp = []
	#for bins_dat,bins_tim in zip(res_bins_dat,res_bins_tim):
	  #dum,tind = filterTimestamps(tbins.time,bins_tim)
	  #res_bins_tim_tmp.extend([bins_tim[ttind] for ttind in tind])
	  #res_bins_dat_tmp.extend([bins_dat[ttind] for ttind in tind])
	#res_bins_dat = res_bins_dat_tmp
	#res_bins_tim = res_bins_tim_tmp
	#res_binvec.append(tbin.scan[tbin.scan.keys()[0]])
      #else:
	#res_bins_dat_tmp = []
	#res_bins_tim_tmp = []
	#for bins_dat,bins_tim in zip(res_bins_dat,res_bins_tim):
	  #inds,bins = digitize(bins_dat,tbins)
	  ## continue here ToDo
	  #binrange = range(1,len(bins))
	  #tind = [(inds==tbin).nonzero()[0] for tbin in binrange]
	  #res_bins_tim_tmp.extend([bins_tim[ttind] for ttind in tind])
	  #res_bins_dat_tmp.extend([bins_dat[ttind] for ttind in tind])
	#res_bins_dat = res_bins_dat_tmp
	#res_bins_tim = res_bins_tim_tmp
      #return memdata(input=[res_bins_dat,res_bins_tim])


    #####


    if type(bins) is tuple:
      Ndim = len(bins)
      for tbins in bins:
	if isinstance(tbins,memdata):
	  tb = tbins.ones()*np.arange(len(tbins))
	  tbInd,stsz = ravelScanSteps(tb.data)
	  tbTim,stsz = ravelScanSteps(tb.time)



  else:
    if not figName:
      figName = 'Select filter limits'

    if graphicalInput and bins==None:
      tools.nfigure(figName)
      pl.clf()
      N,edg = tools.histogramSmart(dat)
      pl.step(edg[:-1]+np.diff(edg),N,'k')
      lims = tools.getSpanCoordinates()
      ip = 0
      bins = None
      hs = []
      while not ip=='q':
	ip = raw_input('Enter number of bins (q to finish)')
	if ip=='q': continue
	bins = np.linspace(np.min(lims),np.max(lims),int(ip))
	for th in hs: 
	  th.remove()
	  hs = []
	hs = [pl.axvline(te) for te in bins]
      print "Selected bins: %g, %g,  "%(np.min(lims),np.max(lims))


  return np.digitize(dat,bins),bins



    

 



def parameterFilt(par,dataset=None,name=None,lims=None,graphicalInput=True,scanstep=None,figName=None):
  par = iterfy(par)
  if dataset:
    dataset.filtTimestamps()

  if not scanstep:
    dat = np.hstack([spar for spar in par])
  else:
    dat = par[scanstep]
  if not figName:
    figName = 'Select filter limits'

  if graphicalInput and not lims:
    tools.nfigure(figName)
    pl.clf()
    N,edg = tools.histogramSmart(dat)
    pl.step(edg[:-1]+np.diff(edg),N,'k')
    lims = tools.getSpanCoordinates()

  tfilt = []
  for tpar in par:
    tfilt.append(~((tpar>lims[0])&(tpar<lims[1])))
  if dataset:
    #if not name:
      #name = 'unnnamed'
    #filtname = '_filt_'+name
    #dataset.__dict__[filtname] = tfilt
    #dataset._initfilter()
    #return lims
    dataset.filter = tfilt
    dataset._initfilter()
    return lims
  else:
    return tfilt,lims

def corrFilt(par0,par1,lims=None,graphicalInput=True,scanstep=None,figName=None,ratio=False):
  
 

  if not figName:
    figName = 'Select filter limits'

  if graphicalInput and not lims:
    tools.nfigure(figName)
    pl.clf()
    if ratio:
      pl.plot(par0.R,par1.R/par0.R,'.k',ms=1)
    else:
      pl.plot(par0.R,par1.R,'.k',ms=1)
    lims = tools.getRectangleCoordinates()
    lims = list(np.reshape(lims,[2,-1]))
    print 'chosen limits are %s' %lims

  filt0 = par0.filter(lims=list(lims[0])).ones()
  filt1 = par1.filter(lims=list(lims[1])).ones()

  filt = filt0*filt1


  return filt,lims


def getRunInfo(run,epicsPVs=[]):
  try:
    d = dataset(run)
    output = dict()
    output['detectors'] = d.detectors
    time = d.__dict__[d.detectors[0]].time
    output['noofsteps'] = len(time)
    output['shots_per_step'] = [len(time[n]) for n in range(output['noofsteps'])]
    output['starttime'] = datetime.datetime.fromtimestamp(time[0][0][0])
    output['endtime'] = datetime.datetime.fromtimestamp(time[-1][-1][0])
    if 'scanMot' in d.__dict__.keys():
      output['scanmot'] = d.__dict__['scanMot']
    else:
      output['scanmot'] = ''

    if 'scanVec' in d.__dict__.keys():
      output['scanvec'] = d.__dict__['scanVec']
    else:
      output['scanvec'] = []

    pvs = []
    for PV in epicsPVs:
      if 'PVnames' in d.epics.__dict__:
        if PV in d.epics.PVnames:
          pvs.append((PV,d.epics.getEpicsPVdata(PV)[0][0]['value']))
        else:
          pvs.append((PV,[]))
      else:
        pvs.append((PV,[]))
    output['epics'] = pvs


    return output,d

  except:
    print "Could not load run"
    print run

def printRunInfo(exp,runnos,epicsPVs=[],mode='ascii'):
  totlist = []
  for run in runnos:
    try:
      info,d = getRunInfo((exp,run),epicsPVs)
      titles   = []
      outplist = []
      titles.append('Run')
      outplist.append('%d'%(run))
      
      titles.append('Start time')
      outplist.append(info['starttime'].time().isoformat())

      titles.append('Steps')
      outplist.append(str(info['noofsteps']))

      titles.append('Events')
      outplist.append('%d'%(np.sum(info['shots_per_step'])))
      
      titles.append('Motor')
      outplist.append(str(info['scanmot']))
      
      svec = info['scanvec']

      titles.append('Scan range')
      titles.append('Stepsize')
      if not svec==[]:
        outplist.append('%g to %g'%(np.min(svec),np.max(svec)))
        if np.mean(np.diff(svec,n=2))/np.mean(np.diff(svec,n=1)) < .001:
          outplist.append('%g'%(np.mean(np.diff(svec))))
        else:
          outplist.append('')
      else:
        outplist.append('')
        outplist.append('')


      for pv in info['epics']:
        titles.append(pv[0])
        if not pv[1] ==[]:
          outplist.append('%g'%(pv[1]))
        else:
          outplist.append('')

      titles.append('Detectors')
      outplist.append(', '.join(info['detectors']))
      totlist.append(outplist)
    except:
      print 'Failure with Run %d'%(run)

  if mode is 'ascii' or 'print':
    maxsz = []
    itlist = list(totlist)
    itlist.append(titles)
    for nfield in range(len(totlist[0])):
      maxsz.append(np.max([len(outplist[nfield]) for outplist in itlist]))
    fmt = ''
    for maxchar in maxsz:
      fmt += '%'+str(maxchar)+'s|'
    outpstring = ''
    outpstring += fmt %tuple(titles)
    outpstring += '\n'
    for outplist in totlist:
      outpstring += fmt %tuple(outplist)
      outpstring += '\n'
    print outpstring

  return titles,outplist,outpstring

def saveRunsToCache(exp,runnos,dets=[]):
  for runno in runnos:
    try:
      d = dataset((exp,runno),dets)
      d.save(force=True)
    except:
      print "Could not save run %d in %s" %(runno,exp)




  

def ipmfilt(ipmdet):
  d = [np.asarray([c0,c1,c2,c3]).transpose() for c0,c1,c2,c3 in zip(ipmdet._channel0,ipmdet._channel1,ipmdet._channel2,ipmdet._channel3)]
  dall = np.vstack(d)
  ph = []
  f = tools.nfigure('ipmfilt')
  pl.clf()
  hy0,hx0 = np.histogram(dall[:,0],bins = np.linspace(np.log10(min(dall[:,0])),np.log10(max(dall[:,0])),len(dall[:,0])/30))
  #de=bug
  hy2,hx2 = np.histogram(dall[:,2],bins = np.logspace(np.log10(min(dall[:,2])),np.log10(max(dall[:,2])),len(dall[:,0])/30))
  ph.append(f.add_axes([.1,.4,.35,.5]))
  ph[-1].plot(dall[:,0],dall[:,1]/dall[:,0],'.k',ms=1)
  ph[0].set_xscale('log')
  ph[0].set_yscale('log')
  #ph.append(pl.subplot(2,2,2))
  ph.append(f.add_axes([.55,.4,.35,.5]))
  ph[-1].plot(dall[:,2],dall[:,3]/dall[:,2],'.k',ms=1)
  ph[1].set_xscale('log')
  ph[1].set_yscale('log')
  #ph.append(pl.subplot(2,2,3))
  ph.append(f.add_axes([.1,.1,.35,.3],sharex=ph[0]))
  ph[-1].plot(np.diff(hx0),hy0,'k')
  ph.append(f.add_axes([.1,.1,.35,.3]))
  #ph.append(pl.subplot(2,2,4))
  ph.append(f.add_axes([.55,.1,.35,.3],sharex=ph[1]))
  ph[-1].plot(np.diff(hx2),hy2,'k')


def _rdDetStepData(fina,step,shots,det=''):
  d = dataset(fina,[det])
  data = d.__dict__[det].rdStepData(step,shots)
  return data

def _rdDetAllData(fina,det=''):
  d = dataset(fina,[det])
  data = d.__dict__[det].rdAllData()
  return data

def _compress(*args):
  #print kwargs
  #print args
  #data = kwargs['data']
  data = args[0]
  dout = [dat.compressed() for dat in data]
  return dout

def _getInputDirectory(cnfFile):
    knownhosts = cnfFile['defaultPath'].keys()
    if gethostname() in knownhosts:
      inputDirectory = cnfFile['defaultPath'][gethostname()]
    else:
      inputDirectory = os.path.join(cnfFile['defaultPath']['default'],cnfFile['beamline'])
    return inputDirectory



def getValidFieldname(name,lowerit=True):
  if lowerit:
    name = name.lower()
  name = name.replace(':','_')
  name = name.replace('.','_')
  name = name.replace(' ','_')
  name = name.replace('-','_')
  #if lowerit:
    #while name.find('_')>-1:
      #s,e = name.split('_')
      #name = s+e[0].upper()+e[1:]
  return name


class _Lcls_beamline(object):
  def __init__(self,cnfFile):
    self.cnfFile = cnfFile
  def __getattr__(self,inp):
    self.__dir__()
    return self.__getattribute(inp)
  def __dir__(self):
    ipflist = os.listdir(_getInputDirectory(self.cnfFile))
    explist  = []
    explist = [file  for file in ipflist if re.match(self.cnfFile['beamline'] + '[micom0123456789]',file)]
    #for exp in getExperimentList(printit=False)['no']:
      #if exp in ipflist:
        #explist.append(exp)
    self.explist = explist
    output =  copy(explist)
    output.extend(self.__dict__.keys())
    for exp in explist:
      self.__dict__[exp] = _Experiment(exp,self.cnfFile)
    #print output
    return output


class _Experiment(object):
  def __init__(self,exp,cnfFile):
    self.cnfFile = cnfFile
    self.exp = exp
  def __getattr__(self,inp):
    self.__dir__()
    return self.__getattribute(inp)
  def __dir__(self):
    filelist = os.listdir(_getInputDirectory(self.cnfFile) + '/' + self.exp + '/hdf5')
    runfilelist = [file  for file in filelist if re.match(self.exp + '-r[0123456789]*.h5',file)]
    runlistnum = [int(file.split(self.exp + '-r')[1].split('.h5')[0]) for file in runfilelist]
    runlist = ['run%04d'%runno for runno in runlistnum]
    self.runlist = runlistnum
    output =  copy(runlist)
    output.extend(self.__dict__.keys())
    for run in runlistnum:
      runname = 'run%04d'%run
      self.__dict__[runname] = _Run(self.exp,run,beamline=self.cnfFile['beamline'])
    return output
    
class _Run(object):
  def __init__(self,exp,run,beamline=None):
    self.exp = exp
    self.run = run
    self.beamline = beamline
  def __getattr__(self,inp):
    self.__dir__()
    return self.__getattribute(inp)
  def __dir__(self):
    exp = self.exp
    run = self.run
    print "Loading dataset information..."
    d = dataset((exp,run),beamline=self.beamline)
    print "done!"
    self.dataset = d
    detectors = d.detectors
    self.detectors = copy(detectors)
    if 'epics' in d.__dict__.keys():
      detectors.append('epics')
    if 'scanMot' in d.__dict__.keys():
      detectors.append('scanMot')
    if 'scanVec' in d.__dict__.keys():
      detectors.append('scanVec')
    
    output =  copy(detectors)
    output.extend(self.__dict__.keys())
    for det in detectors:
      self.__dict__[det] = d.__dict__[det]
    #print output
    return output

  
############## Dataset finder ###############

def findDatasetsInHdf5(filehandle):
  rundset = filehandle['Configure:0000']['Run:0000']
  ccname = rundset.keys()[0]
  dsets = crawlforDatasets(rundset[ccname])
  return dsets

def parseToCnf(fileHandle):
  cnf = dict()
  cnf['areaDet'] = dict()
  cnf['pointDet'] = dict()
  dsets = findDatasetsInHdf5(fileHandle)

  for dset in dsets:
    detname = tools.varName(os.path.split(os.path.split(dset['dset_data'])[0])[1])
    ccn = '/Configure:0000/Run:0000/CalibCycle:0000/'

    ds = fileHandle[dset['dset_data']]
    if len(np.shape(ds[0]))>0:
      cnf['areaDet'][detname] = dict()
      cnf['areaDet'][detname]['conf'] = 'dummy'
      cnf['areaDet'][detname]['timestamp'] = dset['dset_time'].split(ccn)[1]
      cnf['areaDet'][detname]['data'] = dset['dset_data'].split(ccn)[1]
    else:
      cnf['pointDet'][detname] = dict()
      cnf['pointDet'][detname]['conf'] = 'dummy'
      cnf['pointDet'][detname]['timestap'] = dset['dset_time'].split(ccn)[1]
      cnf['pointDet'][detname]['data'] = dset['dset_data'].split(ccn)[1]
  return cnf
  
def crawlforDatasets(group,skipEpics = True, skipEvr = True):
  found = []
  if 'Epics::EpicsPv' in group.name:
    return found
  if 'Evr' in group.name:
    return found
  itemnames  = group.keys()
  if u'time' in itemnames:
    if isinstance(group['time'],h5py.Dataset):
      if u'image' in itemnames:
        if isinstance(group['image'],h5py.Dataset):
          if group['time'].shape[0] == group['image'].shape[0]:
            dset = dict()
            dset['dset_data'] =  group['image'].name
            dset['dset_time'] =  group['time'].name
            #dset['dset_type'] =  'areadet'
            found.extend([dset])
      elif u'data' in itemnames:
        if isinstance(group['data'],h5py.Dataset):
          if group['time'].shape[0] == group['data'].shape[0]:
            dset = dict()
            dset['dset_data'] =  group['data'].name
            dset['dset_time'] =  group['time'].name
            #dset['dset_type'] =  'pointdet'
            found.extend([dset])
  else:
    for itemname in itemnames:
      if isinstance(group[itemname],h5py.Group):
        found.extend(crawlforDatasets(group[itemname]))
  return found

###### WRAPPING #######

corrNonlin = wrapFunc(tools.corrNonlin)

############## GENERAL STUFF ###############

beamlines = ['amo','sxr','xpp','xcs','cxi','mec','mob']

cnfFiles = dict()
for bl in beamlines:
  cnfFiles[bl] = rdConfiguration(beamline=bl)
  exec(bl+'=_Lcls_beamline(cnfFiles[bl])')

cnfFile = rdConfiguration()
beamline = cnfFile['beamline']
# initialize pointdet readers
#for det in cnfFile['pointDet'].keys():
#  exec('rd'+det[0].upper()+det[1:]+'AllData = partial(_rdDetAllData,det=det)')

# initialize areadet readers
for det in cnfFile['areaDet'].keys():
  exec('rd'+det[0].upper()+det[1:]+'StepData = partial(_rdDetStepData,det=det)')

point_detectors = cnfFile['pointDet'].keys()
area_detectors = cnfFile['areaDet'].keys()
detectors = point_detectors
detectors.extend(area_detectors)


