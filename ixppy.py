import pylab as pl
import h5py
import os
from os.path import expanduser
from socket import gethostname
import dateutil
import sys
import numpy as np
import ixppy_specialdet
import tools
from functools import partial
from copy import copy
import re
import examples
import ixppy
import datetime
import operator
#import cspad

# Config file format
# Default instrument
# Point detectors:
# Alias, dataset-data, dataset-timestamp, dataset-configuration 1, dataset-configuration 2, isIpmOrSimilarQuaddiode, hasFexCorrection
# Area detectors
# Alias, dataset-data, dataset-timestamp, dataset-configuration 1, dataset-configuration 2, slabsize, correctionmethods
# Directories
# default data directory (like /red/d/psdm/), main output directory (scratch/or local folder to "take data home").


#print os.path.dirname(__file__)

#print os.path.abspath(__file__)
############## DATASET ############################

class dataset(object):
  def __init__(self,
    inputFileOrExpRunTuple='',
    detectors = [],
    beamline = None,               
    rdPointDetectorsImmediately=True,
    rdTimestamps=True,
    sortForTimestamps=True,
    inputDirectory = 'config',
    outputDirectory = '',
    readCachedData = True,
    cacheFile = None,
    ignore_daqfile = False
              ):
    # init start
    self.config = dropObject()
    if not beamline==None:
      self.config.beamline = beamline
      self.config.cnfFile = rdConfiguration(beamline=beamline)
    else:
      self.config.cnfFile = rdConfiguration()
      self.config.beamline = self.config.cnfFile['beamline']
    self.config.hostname = gethostname()
    self.config.pointDetectors = []
    self.config.areaDetectors = []
    #self._cache_directory = 
    self._getInputDirectory(inputDirectory)
    self._getFilename(inputFileOrExpRunTuple)
    if ignore_daqfile:
      self.config.daqfile_available = False
    self.config.fileHandle = []
    self.pointDet = dropObject()
    self.areaDet  = dropObject()
    self.config.base = self
    self.config.scanDataSets = []
    self._controlPv = []
    self.noofsteps = []
    self._filter = []
    self._saved_objects = []
    self._savelist = ['_savelist','_saved_objects']
    #self._dividePointArrayDets()
    #self._initDetectorReaders()
    if self.config.daqfile_available:
      if detectors==[]:
        detectors = self._find_detectors()
      elif detectors=='parse':
        if not self.config.fileHandle:
          self.config.fileHandle = h5py.File(self.config.filename[0],mode='r')
        print "parsing detectors in hd5 file"
        cnf = parseToCnf(self.config.fileHandle)
        self.config.cnfFile = tools.dict_merge(self.config.cnfFile,cnf)
        detectors = cnf['areaDet'].keys()+cnf['pointDet'].keys()

    if not outputDirectory:
      if self.config.cnfFile['cache_directory']:
        self.config.outputDirectory = self.config.cnfFile['cache_directory'][0]
      else:
        self.config.outputDirectory = os.getcwd()
    else:
      self.config.outputDirectory = outputDirectory

    if not cacheFile:
      tfina = os.path.split(self.config.filename[0])[-1]
      tfina = os.path.join(self.config.outputDirectory,tfina) 
      if os.path.isfile(tfina):
        if not tfina==self.config.filename:
          self.config.cache_filename = tfina
      else:
        self.config.cache_filename = ''
    else:
        self.config.cache_filename = cacheFile

    self.config.cache_fileHandle = []
    self.config.cache = cache(self.config)

    self.detectors = iterfy(detectors)
    self.initdetectors()

    # initialize epics pvs
    self.epics = epics(self.config)
    self._saved_objects.append('epics')
    #print readCachedData
    #print os.path.isfile(self.config.cache_filename)
    if readCachedData and os.path.isfile(self.config.cache_filename):
      self.config.cache.load()
      self._unwrap_detectors()
      self._initfilter()

    if rdPointDetectorsImmediately:
      self.rdPointDetectors()
    if rdTimestamps:
      self.filtTimestamps()
    ##### END INIT #####

  def initdetectors(self):
    for det in self.detectors:
      self._initdet(det)

  def __getitem__(self,x):
    return eval("self.%s"%x)

  def rdPointDetectors(self):
    """reads all point detectors from the defined list of detectors"""
    for det in self.config.pointDetectors:
      data = self.__dict__[det].rdAllData()
      self.__dict__[det]._add_saved_datafield('data',data)
    self._unwrap_detectors()
 
  def _find_detectors(self):
    if not self.config.fileHandle:
      self.config.fileHandle = h5py.File(self.config.filename[0],mode='r')
    dets = []
    for det in self.config.cnfFile['areaDet'].keys():
      odet = det
      det = det.split('_bak')[0]
      if det in dets: continue
      dset_time = self.config.cnfFile['areaDet'][odet]['dataset_time']
      dsetstring = self._mkScanStepDataSetString(dset_time,0)
      try:
        self.config.fileHandle[dsetstring]
        dets.append(det)
      except:
        pass
    for det in self.config.cnfFile['pointDet'].keys():
      odet = det
      det = det.split('_bak')[0]
      if det in dets: continue
      dset_time = self.config.cnfFile['pointDet'][odet]['dataset_time']
      dsetstring = self._mkScanStepDataSetString(dset_time,0)
      try:
        self.config.fileHandle[dsetstring]
        dets.append(det)
      except:
        pass
    return dets


  def _unwrap_detectors(self):
    """get structured array components from detector data and mke them masked arrays for easier use"""
    for det in self.config.pointDetectors:
      self.__dict__[det]._unwrap_data()

  def rdTimestamps(self):
    """reads timestamps of all selected detectors"""
    for det in self.detectors:
      self.__dict__[det]._add_saved_datafield('time',self._rdDetTimestamp(det))

  def filtTimestamps(self):
    """tests the detector timestamps for equality and generates a filter for each selected detector"""
    # check if timestamps are already read.
    for det in self.detectors:
      if not hasattr(self.__dict__[det],'time'):
	print "reading timestamps of %s" %det
        self.rdTimestamps()
        break
      
    for det in self.detectors:
      self.__dict__[det]._add_saved_datafield('_filt_time',[])
    for sNo in range(self.noofsteps):
      for det in self.detectors:
        try:
          tfilt = np.ones(np.shape(self.__dict__[det].time[sNo]),dtype=bool)
          ttime = self._timeInMs(self.__dict__[det].time[sNo])
          for cdet in np.setdiff1d(self.detectors,[det]):
            ctime = self._timeInMs(self.__dict__[cdet].time[sNo])
            cfilt = np.in1d(ttime,ctime)
            #if sum(np.logical_not(cfilt))>0:
              #de=bug
            tfilt = np.logical_and(tfilt,cfilt)
          tfilt = np.logical_not(tfilt)
          self.__dict__[det]._filt_time.append(tfilt)
        except:
          print "Problem in time stamp filtering with %s at step %d." %(det,sNo)
    self._initfilter()

  def old_initfilter(self):
    for det in self.detectors:
      # get all filter names
      timefilt = []
      for k in self.__dict__[det].__dict__.keys():
        if '_filt_time' in k:
          timefilt.append(k)
      # make global filter
      filt_tot = []
      for filter in filters:
        if not filt_tot:
          filt_tot = self.__dict__[det].__dict__[filter]
        for sNo in range(len(filt_tot)):
          filt_tot[sNo] = np.logical_and(filt_tot[sNo],self.__dict__[det].__dict__[filter][sNo])
      # add total filter to masked arrays
      for sNo in range(len(filt_tot)):
        for name in self.__dict__[det].data[sNo].dtype.names:
          self.__dict__[det].__dict__[name][sNo].mask = np.logical_not(filt_tot[sNo])

  #def _initfilter(self,filtername=None):
    #for det in self.detectors:
      ## get all filter names, timefilter first, then global filters
      #allfilters = []
      #for k in self.__dict__[det].__dict__.keys():
        #if '_filt_time' in k:
          #allfilters.append(self.__dict__[det].__dict__[k])
      ##print allfilters
      
      #if not filtername:
        ## Get global filters, first names, then names are sorted then the filter data. could be sorted by size at some point.
        #allGfilters = []
        #for k in self.__dict__.keys():
          #if '_filt_' in k:
            #allGfilters.append(k)
        #allGfilters.sort()
        #for key in allGfilters:
          #allfilters.append(self.__dict__[key])
        ##print allfilters
      #else:
        #filtername = iterfy(filtername)
        #allGfilters = [self.__dict__['_filt_'+tfn] for tfn in filtername]
      
      #filt_tot = []
      #for filter in allfilters:
        #if not filt_tot:
          #for sNo in range(self.noofsteps):
            #filt_tot.append(filter[sNo].copy())
        #else:
          #for sNo in range(len(filt_tot)):
            #tunmasked = filt_tot[sNo][~filt_tot[sNo]].copy()
            #tunmasked[filter[sNo]] = True
            #filt_tot[sNo][~filt_tot[sNo]] = tunmasked.copy()
      #self.__dict__[det].filter = filt_tot
      #try:
        ## add total filter to masked arrays
        #for sNo in range(len(filt_tot)):
          ##for name in self.__dict__[det].data[sNo].dtype.names:
          #for name in self.__dict__[det]._masked_arrays['data']:
            ##print len(filt_tot[sNo])
            #self.__dict__[det].__dict__[name][sNo].mask = filt_tot[sNo]
      #except Exception,e:
        ##print e
        ##print "problem in initfilter"
        #pass
  ######quick hack for global filter
    #allfilters = []
    #for key in allGfilters:
      #allfilters.append(self.__dict__[key])
    #filt_tot = []
    #for filter in allfilters:
      #if not filt_tot:
        #for sNo in range(self.noofsteps):
          #filt_tot.append(filter[sNo].copy())
      #else:
        #for sNo in range(len(filt_tot)):
          #tunmasked = filt_tot[sNo][~filt_tot[sNo]].copy()
          #tunmasked[filter[sNo]] = True
          #filt_tot[sNo][~filt_tot[sNo]] = tunmasked.copy()
    #self.filter = filt_tot

  def _mergefilters(self):
    filter = []
    for sNo in range(self.noofsteps):
      tstepf = np.vstack([tfilt[sNo] for tfilt in self._filter])
      filter.append(tstepf.any(axis=0))
    return filter


  def _addfilter(self,newfilter):
    mergedfilter = self.filter
    tfilter = self._filter[0]
    Ofilter = []
    for Stfilter,Snewfilter,Smergedfilter in zip(tfilter,newfilter,mergedfilter):
      if len(Stfilter)==len(Snewfilter):
        Ofilter.append(Snewfilter)
      elif sum(~Smergedfilter)==len(Snewfilter):
        tmpf = Smergedfilter
        tunmasked = tmpf[~tmpf].copy()
        tunmasked[Snewfilter] = True
        tmpf[~tmpf] = tunmasked.copy()
        Ofilter.append(tmpf)
      else:
        print "NB: Filter step doesn't fit in size!"
    self._filter.append(Ofilter)
    self._initfilter()

  filter = property(_mergefilters,_addfilter)

  def filterClear(self,remove=-1):
    if type(remove) is str and remove=='all':
      self._filter = []
    else:
      remove = iterfy(remove)
      for ind in remove:
	self._filter.pop(ind)
    self._initfilter()

      

  def _initfilter(self,filtername=None):

    if not self._filter:
      
      # get all detector filter, check if all are in line
      alldetfilters = []
      for det in self.detectors:
        for k in self.__dict__[det].__dict__.keys():
          if '_filt_time' in k:
            alldetfilters.append(self.__dict__[det].__dict__[k])

      tfilter = []
      for sNo in range(self.noofsteps):
        lengths = np.array([np.sum(~detfilt[sNo]) for detfilt in alldetfilters])
        if ~(np.sum(np.abs(np.diff(lengths)))==0):
          print "something went wrong with timestampfiltering in scan step %d."%(sNo)
          #de=bug
        tfilter.append(np.zeros(lengths[0],dtype=bool))
      self._filter.append(tfilter)
    
    # Now initialize filter for each detector 
    for det in self.detectors:
      # get all filter names, timefilter first, then global filters
      for k in self.__dict__[det].__dict__.keys():
        if '_filt_time' in k:
          detfilter = self.__dict__[det].__dict__[k]

      Gfilter = self.filter
      filt_tot = [tdf.copy() for tdf in detfilter]
      for sNo in range(len(filt_tot)):
        tunmasked = filt_tot[sNo][~filt_tot[sNo]].copy()
        tunmasked[Gfilter[sNo]] = True
        filt_tot[sNo][~filt_tot[sNo]] = tunmasked.copy()
      self.__dict__[det].filter = [ft.copy() for ft in filt_tot]
      # add total filter to masked arrays
      for sNo in range(len(filt_tot)):
        #for name in self.__dict__[det].data[sNo].dtype.names:
        for name in self.__dict__[det]._masked_arrays['data']:
          #print len(filt_tot[sNo])
          self.__dict__[det].__dict__[name][sNo].mask = filt_tot[sNo]


  def _add_saved_datafield(self,name,data):
    """adds to a list of datasets to be saved. For faster data recovery and for saving custom analysis progresses (e.g. reduced data from pixel detectors)."""
    if not "_savelist" in self.__dict__:
      self._savelist = ['_savelist']
    self._savelist.append(name)
    self._savelist = list(set(self._savelist))
    self.__dict__[name] = data


  def _timeInMs(self,time):
    """Makes millisecond array from the read second and nanosecond arrays"""
    ms = np.uint64(time['seconds'])
    ms = ms*1000 + np.round_(time['nanoseconds'],-6)/1000000
    return ms

  def _initdet(self,det):
    class tclass(singledet):
      pass

    tdet = tclass(self.config,det)
    self.__dict__[det] = tdet
    self._saved_objects.append(det)

    if self.__dict__[det]._ispointdet():
      self.pointDet.__dict__[det] = tdet
      self.config.pointDetectors.append(det)
    if self.__dict__[det]._isareadet():
      self.areaDet.__dict__[det] = tdet
      self.config.areaDetectors.append(det)

  def _rdScanPar(self):
    """Reads the scan datasets and the scanned motor(s) into the dataset structure."""
    if not self.config.fileHandle:
      self.config.fileHandle = h5py.File(self.config.filename[0],mode='r')
    if not self.config.scanDataSets:  
      self.config.scanDataSets = self.config.fileHandle[self.config.cnfFile['scan_step'][0]].keys()
    if not self.noofsteps:  
      self._add_saved_datafield('noofsteps',len(self.config.scanDataSets))
    if not self._controlPv:
      controlPv = []
      for tdsetbas in self.config.scanDataSets:
        tdsetnames = [os.path.join(self.config.cnfFile['scan_step'][0],tdsetbas,'ControlData::ConfigV1/Control/pvControls'),
                      os.path.join(self.config.cnfFile['scan_step'][0],tdsetbas,'ControlData::ConfigV1/NoDetector.0/pvControls'),
                      os.path.join(self.config.cnfFile['scan_step'][0],tdsetbas,'ControlData::ConfigV2/Control/pvControls')]
        trycycle = 0
        done=False
        while trycycle<len(tdsetnames) and not done:
          try:
            tdsetname = tdsetnames[trycycle]
            try:
              controlPv.append(self.config.fileHandle[tdsetname].value)
            except:
              controlPv.append(self.config.fileHandle[tdsetname])
            done = True
          except:
            trycycle+=1
            continue
      self._controlPv = controlPv
      try:
        self._getScanMotVec()
      except:
        pass
      return

  def _getScanMotVec(self):
    names = self._controlPv[0]['name']
    if len(names)>1:
      tname = [name for name in names]
      self._add_saved_datafield('scanMot',tname)
      scanVec=[]
      for motNo in range(len(tname)):
        tscanVec = []
        for val in self._controlPv:
          tscanVec.append(val['value'][motNo])
        scanVec.append(np.array(tscanVec))
      self._add_saved_datafield('scanVec',scanVec)
    else:
      self._add_saved_datafield('scanVec',np.array([val['value'][0] for val in self._controlPv]))
      tname = self._controlPv[0]['name'][0]
      #print tname
      self._add_saved_datafield('scanMot',tname)

  def _mkScanStepDataSetString(self,dsetString,stepIndex):
    """Makes dataset sting for certain scan step and certain dtector daset (the fraction of the path string after the "calib cycle")"""
    string0 = os.path.join(self.config.cnfFile['scan_step'][0],\
             self.config.cnfFile['scan_step'][1])
    #digits  = int(self.config.cnfFile['scan_step'][2]) 
    output = []
    for ind in iterfy(stepIndex):
      ind = '%04d' %ind
      output.append(os.path.join(string0+ind,dsetString))
    if len(output) is 1:
      output = output[0]
    return output
      

    

  def _rdDetTimestamp(self,detector):
    self._rdScanPar()
    for tdetdset in self.__dict__[detector]._dataset_time:
      try:
        data = ixppyList()
        for stepNo in range(len(self._controlPv)):
          dsetstring = self._mkScanStepDataSetString(tdetdset,stepNo)
          #print dsetstring
          dset = self.config.fileHandle[dsetstring]
          data.append(rdHdf5dataFromDataSet(dset))
        #print "found time data for %s" % detector
        return data
        #break
      except Exception,e:
        #print e
        continue
        #try:
          #return data
          #continue
        #except:
          #continue
      
  def _rdDetAllData(self,detector):
    self._rdScanPar()
    for tdetdset in self.__dict__[detector]._dataset_data:
      try:
        data = ixppyList()
        for stepNo in range(len(self._controlPv)):
          dsetstring = self._mkScanStepDataSetString(tdetdset,stepNo)
          dset = self.config.fileHandle[dsetstring]
          data.append(rdHdf5dataFromDataSet(dset))
        return data
        break
      except:
        #print "Problem reading %s at step %d !!!" %(detector,stepNo)
        continue
        #try:
          #return data
        #except:
          #print "Problem reading %s at step %d !!!" %(detector,stepNo)
          #continue
        #except:
          #continue


  def _rdDetStepData(self,detector,stepNo=0,shotNos='all'):
    self._rdScanPar()
    #tdetdset = self.__dict__[detector]._dataset_data
    for tdetdset in self.__dict__[detector]._dataset_data:
      try:
        dsetstring = self._mkScanStepDataSetString(tdetdset,stepNo)
        dset = self.config.fileHandle[dsetstring]
        data = rdHdf5dataFromDataSet(dset,shotNos)
        return data
        break
      except:
        continue

    
  def all_detectors(self):
    """Prints list of all possible detectors in the present configuration. Useful when not certain about the aliases to use."""
    pointdets = ''
    for i in self.config.cnfFile['pointDet'].keys():
      if not '_bak' in i:
        pointdets+= ', ' + i
    pointdets = pointdets.strip(', ')
    print "Point detectors\n %s" %(pointdets)
    print '\n'
    areadets = ''
    for i in self.config.cnfFile['areaDet'].keys():
      if not '_bak' in i:
        areadets+= ', ' + i
    areadets = areadets.strip(', ')
    print "Area detectors\n %s" %(areadets)



  def _getInputDirectory(self,ipd):
    if ipd=='config':
      knownhosts = self.config.cnfFile['defaultPath'].keys()
      if self.config.hostname in knownhosts:
        self.config.inputDirectory = self.config.cnfFile['defaultPath'][self.config.hostname]
      else:
        self.config.inputDirectory = os.path.join(self.config.cnfFile['defaultPath']['default'],self.config.beamline)

  def _getFilename(self,input):
    if type(input) is str:
      filenames = input
    if type(input) is tuple:
      if type(input[0]) is str:
        self.config.experiment = input[0]
      if type(input[0]) is int:
        self.config.experiment = self.config.beamline+str(input[0])
      self.config.run = input[1]
      filenames = self._makeFilenameFromExpRun()
    self.config.filename = iterfy(filenames)
    for f in self.config.filename:
      if (not os.path.exists(f)):
        print "Asked to read file %s, but is does not exist" % f
	self.config.daqfile_available = False
      else:
	self.config.daqfile_available = True
       

  def _makeFilenameFromExpRun(self):
    if type(self.config.run) is not list:
      run = [self.config.run]
    else:
      run = self.config.run
    filenames = []
    for trun in run:
      tpath = '%s/hdf5'  %(self.config.experiment)
      tfile = '%s-r%04d.h5'  %(self.config.experiment,trun)
      filenames.append(os.path.join(self.config.inputDirectory,tpath,tfile))
    return filenames


  def _rdCachedData(self):
    fina = self.config.cache_filename
    F = h5py.File(fina)
    fkeys = F.keys()

  def addSavedObject(self,name):
    self.__dict__[name] = dropData()
    self.__dict__[name]._savelist = []
    self._saved_objects.append(name)
    
  def _add_saved_datafield(self,name,data):
    if not "_savelist" in self.__dict__:
      self._savelist = ['_savelist']
    self._savelist.append(name)
    self._savelist = list(set(self._savelist))
    self.__dict__[name] = data

  def save(self,force=False):
    self.config.cache.save(force=force)
      
############## DETECTOR ############################

class singledet(object):
  def __init__(self,config,det):
    det = det.split('_bak')[0]
    self._detector = det
    for detkind in ['pointDet','areaDet']:
      if det in config.cnfFile[detkind]:
        self._detClass = detkind
        break
    dets = []
    for i in config.cnfFile[self._detClass].keys():
      if det==i.split('_bak')[0]: dets.append(i)

    self._dataset_data = [] 
    self._dataset_time = []
    self._dataset_conf = []
    
    self._masked_arrays = dict()
    self._masked_arrays['data'] = []


    #self._savelist = ['_savelist','_detector','_detClass','_masked_arrays']
    self._savelist = ['_savelist','_detector','_detClass']

    for det in dets:
      self._dataset_data.append(config.cnfFile[self._detClass][det]['dataset_data'])
      self._dataset_time.append(config.cnfFile[self._detClass][det]['dataset_time'])
      self._dataset_conf.append(config.cnfFile[self._detClass][det]['dataset_conf'])
    self.config = config
    if 'special' in self._dataset_data:
      self._specialdet = ixppy_specialdet.__dict__[det.split('_bak')[0]](config)
    else:
      self._specialdet = []

    #if self._detClass = 'areaDet':


  #def _get_h5_dsets(self):
    #for tdset in self._dataset_data:
      #dsetstring = self._mkScanStepDataSetString(tdetdset,0)
      #self.data = self.config.fileHandle[self._dataset
        #dset = self.config.fileHandle[dsetstring]
        #data = rdHdf5dataFromDataSet(dset,shotNos)
        #return data




  def _add_saved_datafield(self,name,data):
    if not "_savelist" in self.__dict__:
      self._savelist = ['_savelist']
    self._savelist.append(name)
    self._savelist = list(set(self._savelist))
    self.__dict__[name] = data

  def append_h5_datalist(self,name,data):
    """Adds data to a datafield (in the detector object)
    that is a list of hdf5 dataset pointers. Each time this 
    function is used a dataset with the defined data is 
    added to the cache file. This should acilitate saving larger 
    datasets with the data without overflowing memory."""
    if not "_h5list" in self.__dict__:
      self._h5list = []
    if not "_savelist" in self.__dict__:
      self._savelist = ['_savelist']
    if not "_h5list" in self._savelist:
      self._savelist.append('_h5list')
    if not name in self._h5list:
      self._h5list.append(name)
    self._h5list = list(set(self._h5list))
    if not hasattr(self,name):
      self.__dict__[name] = []
    oldlen = len(self.__dict__[name])
    dsname = '#'+'%06d' %(oldlen)
    
    # get filehandle
    if not self.config.cache_filename:
      tfina = os.path.split(self.config.filename[0])[-1]
      tfina = os.path.join(self.config.outputDirectory,tfina) 
      self.config.cache_filename = tfina
    fh = self.config.cache.cacheFileHandle

    if self._detector not in fh.keys():
      fh.require_group(self._detector)
    if name not in fh[self._detector].keys():
      gh = fh[self._detector].require_group(name)
    self.__dict__[name].append(
      fh[self._detector][name].create_dataset(dsname,data=data))

    
  
  def rdAllData(self,force=False):
    if hasattr(self,'data') and not force:
      print "%s data exists. Force overwrite with \'force\' option in det.rdAllData." %(self._detector)
      data = self.data
    else:
      if not self._specialdet:
        data  = self.config.base._rdDetAllData(self._detector)
      else:
        data  = self._specialdet.rdAllData()
    return data


  def rdStepData(self,stepNo=0,shotNos='all'):
    if not self._specialdet:
      data  = self.config.base._rdDetStepData(self._detector,stepNo,shotNos)
    else:
      data  = self._specialdet.rdStepData(stepNo,shotNos)
    return data

  def _ispointdet(self):
    if self._detClass is 'pointDet':
      ispoint = True
    else:
      ispoint = False
    return ispoint

  def _isareadet(self):
    if self._detClass is 'areaDet':
      isarea = True
    else:
      isarea = False
    return isarea
  
  def _compress(self,*args):
    #print kwargs
    #print args
    #data = kwargs['data']
    #de=bug
    if len(args) is 1:
      data = args[0]
    elif len(args) is 2:
      data = args[1].data
    dout = [dat.compressed() for dat in data]
    return dout
  
  def _compress_name(self,*args):
    #print kwargs
    #print args
    #data = kwargs['data']
    #de=bug
    if len(args) is 1:
      data = self.__dict__['_'+args[0]]
    elif len(args) is 2:
      data = args[1].__dict__['_'+args[0]]
    dout = [dat.compressed() for dat in data]
    dout = ixppyList(dout)
    return dout

  def _unwrap_data(self):
    """get structured array components from detector data and mke them masked arrays for easier use"""
    data = self.data 
    if data[0].dtype.names:
      for fieldname in data[0].dtype.names:
        if data[0][fieldname].ndim==1:
          self._masked_arrays['data'].append('_'+fieldname)
          self.__dict__['_'+fieldname] = ixppyList() 
          for i in range(len(data)):
            self.__dict__['_'+fieldname].append(np.ma.MaskedArray(data[i][fieldname].copy()))
          self.__dict__['__'+fieldname] = partial(self._compress_name,fieldname)
          setattr(self.__class__,fieldname,property(self.__dict__['__'+fieldname]))
        elif data[0][fieldname].ndim==2:
          noofvecs = np.shape(data[0][fieldname])[1]
          for n in range(noofvecs):
            strfmt = '%0'+'%dd' %(1+np.floor(np.log10(noofvecs)))
            tname = fieldname+strfmt %(n)
            self._masked_arrays['data'].append('_'+tname)
            self.__dict__['_'+tname] = ixppyList()

            for i in range(len(data)):
              self.__dict__['_'+tname].append(np.ma.MaskedArray(data[i][fieldname][:,n].view()))
            self.__dict__['__'+tname] = partial(self._compress_name,tname)
            setattr(self.__class__,tname,property(self.__dict__['__'+tname]))

          
        else:
          print "No clue how to unwrap data %d in %s" %(name,det)

  ########### CHUNK STUFF ###################
  def chunks(self,pointdets=None,Nmax=None,NstepMax=None,steps=None):
    if not pointdets:
      pointdets = dict()
      for det in self.config.pointDetectors:
        #data = self.config.base.__dict__[det].data 

        #if data[0].dtype.names:
          #for name in data[0].dtype.names:
            #tdset = [d for d in self.config.base.__dict__[det].__dict__[name]]
            #pointdets[det + '_' + name] = tdset
        if self.config.base.__dict__[det]._masked_arrays['data']:
          for pdname in self.config.base.__dict__[det]._masked_arrays['data']:
            tdset = self.config.base.__dict__[det].__dict__['_'+pdname]()
            pointdets[det + pdname] = tdset



    self.config.base._rdScanPar()
    filter = self.filter
    indices = []
    for filt in filter:
      indices.append(np.nonzero(~filt)[0])
    
    timestamps = [ts[fli] for ts,fli in zip(self.time,indices)]
    pointdets['time'] = timestamps

    allchunks = []
    if not NstepMax:
      NstepMax = self.config.base.noofsteps
    if not steps:
      steps = range(min(NstepMax,self.config.base.noofsteps))
    for sNO in steps:
      # get dataset
      for tdetdset in self._dataset_data:
        try:
          dsetstring = self.config.base._mkScanStepDataSetString(tdetdset,sNO)
          h5dset = self.config.fileHandle[dsetstring]
          break
        except:
          continue
      if pointdets: 
        tpointdets = dict()
        for tpd in pointdets.keys():
          tpointdets[tpd] = pointdets[tpd][sNO]
      else:
        tpointdets = None
      memavailable = mem().free/8
      allchunks.append(self._chunks(h5dset,indices[sNO],memavailable,Nmax,tpointdets))

    return allchunks
      
  def _chunks(self,dataset,indices,memavailable,maxno,pointdets):
    if not maxno: maxno=np.nan
    elementBytes = dataset.dtype.itemsize
    Nel = min(len(indices),maxno)
    mxsz = np.floor(memavailable/elementBytes)
    Nchunks = int(np.ceil(Nel/mxsz))
    sz = int(np.ceil(1.*Nel/Nchunks))
    chunkind = []
    for n in range(Nchunks):
      chunkind.append(np.arange(n*sz,min((n+1)*sz,Nel)))

    chunks = []
    for cN in range(Nchunks):
      chunks.append(chunk(dataset,indices,chunkind[cN],pointdets))
    return chunks

############## EPICS ############################
class epics(object):
  def __init__(self,config):
    self.config = config
    self._savelist = ['_savelist']
    try:
      self._initEpicsData()
    except:
      print "Epics Data not initialized, possibly not existing in datafile!"
    
    #for PV in self.PVnames:
      #PVflat = PV.replace(':','_')
      #PVflat = PVflat.replace(' ','_')
      #PVflat = PVflat.replace('-','_')
      #PVflat = PVflat.replace('.','_')
      #exec('self.'+PVflat+' = property(self._get_'+PVflat+')')

  def _add_saved_datafield(self,name,data):
    if not "_savelist" in self.__dict__:
      self._savelist = ['_savelist']
    self._savelist.append(name)
    self._savelist = list(set(self._savelist))
    self.__dict__[name] = data

  def _initEpicsData(self):
    self._edset = self.config.cnfFile['epics_dset']
    self._ccedset = self.config.cnfFile['epics_cc_dset']
    if not self.config.fileHandle:
      self.config.fileHandle = h5py.File(self.config.filename[0],mode='r')
    dset = self.config.fileHandle[self._edset]
    pvs = dset.keys()
    self.PVnames = pvs
    for PV in pvs:
      PVflat = PV.replace(':','_')
      PVflat = PVflat.replace('-','_')
      PVflat = PVflat.replace(' ','_')
      PVflat = PVflat.replace('.','_')
      self.__dict__['_get_'+PVflat] = partial(self.getEpicsPVdata,PV)
      setattr(self.__class__,PVflat,property(self.__dict__['_get_'+PVflat]))


  def _rdEpicsPVsingle(self,PV):
    edset = self._edset
    dset = self.config.fileHandle[edset][PV]['data']
    data = dset.value
    return data

    
  def getEpicsPVdata(self,*args):
    #print args
    PV = args[0]      
    PVflat = PV.replace(':','_')
    PVflat = PVflat.replace('-','_')
    PVflat = PVflat.replace(' ','_')
    PVflat = PVflat.replace('.','_')
    PVflat = str(PVflat)

    try:
      data = self.__dict__['_'+PVflat+'_data']
    except:
      data = self._rdEpicsPVdata(PV)
      self._add_saved_datafield('_'+PVflat+'_data',data)
    return data


  def _rdEpicsPVdata(self,PV):
    self.config.base._rdScanPar()
    data = []
    for stepNo in range(len(self.config.base._controlPv)):
      dsetstring = self.config.base._mkScanStepDataSetString(self._ccedset+'/'+PV+'/data',stepNo)
      dset = self.config.fileHandle[dsetstring]
      data.append(rdHdf5dataFromDataSet(dset))
    return data

############ CHUNK #############
class chunk(object):
  def __init__(self,dataset,filtindices,chunkind,pointdets=dict()):
    self.dataset = dataset
    self.filtindices = filtindices
    self.chunkind = chunkind
    if pointdets:
      for pointdet in pointdets.keys():
        self.__dict__[pointdet] = pointdets[pointdet][chunkind]

  def get_data(self):
    shotNos = self.filtindices[pl.ix_(self.chunkind)]
    data = rdHdf5dataFromDataSet(self.dataset,shotNos)
    return data
  data = property(get_data)

  #def get_time(self):
    #shotNos = self.filtindices[pl.ix_(self.chunkind)]
    #data = rdHdf5dataFromDataSet(self.dataset,shotNos)
    #return data
  #time = property(get_time)

############ CACHE #############

class cache(object):
  def __init__(self,config):
    self.config = config
    if os.path.isfile(self.config.cache_filename):
      print 'Found cached data in %s' %(self.config.cache_filename)

  def get_cacheFileHandle(self):
    if not self.config.cache_fileHandle:
      cfh = h5py.File(self.config.cache_filename)
      self.config.cache_fileHandle = cfh
    else:
      cfh = self.config.cache_fileHandle
    return cfh

  def set_cacheFileHandle(self,cfh):
    self.config.cache_fileHandle = cfh

  cacheFileHandle = property(get_cacheFileHandle,set_cacheFileHandle)
 
  def delete(self):
    ristr ='Do you intend to permanently delete file \n  %s\n  (y/n) ' %(self.config.cache_filename)
    if 'y'==raw_input(ristr):
      os.remove(self.config.cache_filename)

  def save(self, cacheFile='default', force=False):
    if cacheFile=='default':
      if not self.config.cache_filename:
        tfina = os.path.split(self.config.filename[0])[-1]
        tfina = os.path.join(self.config.outputDirectory,tfina) 
        self.config.cache_filename = tfina
    else:
      self.config.cache_filename = cacheFile
    
    cfh = self.cacheFileHandle
    # datasets in "root"-dataset instance
    fGroup = cfh.require_group('base_dataset_instance')
    for sfield in self.config.base._savelist:
      if sfield in fGroup.keys() and not force:
        rawstr = 'Overwrite %s in base dataset ? (y/n/a) ' %(sfield)
        ret = raw_input(rawstr)
        if ret=='a': del fGroup[sfield]; force = True
        if ret=='y': del fGroup[sfield]
        if ret=='n': continue
      elif sfield in fGroup.keys() and force:
        print "about to delete %s" %(sfield)
        del fGroup[sfield]
      self.mkDset(fGroup,sfield,self.config.base.__dict__[sfield])

    # detectors and data datasets
    for field in self.config.base._saved_objects:
      fGroup = cfh.require_group(field)
      for sfield in self.config.base.__dict__[field]._savelist:
        if sfield in fGroup.keys() and not force:
          rawstr = 'Overwrite %s in %s ? (y/n/a) ' %(sfield, field)
          ret = raw_input(rawstr)
          if ret=='a': del fGroup[sfield]; force = True
          if ret=='y': del fGroup[sfield]
          if ret=='n': continue
        elif sfield in fGroup.keys() and force:
          print "about to delete %s" %(sfield)
          del fGroup[sfield]
        self.mkDset(fGroup,sfield,self.config.base.__dict__[field].__dict__[sfield])
    cfh.close()
  
  def load(self):
    savobjs = self.cacheFileHandle.keys()
    for savobj in savobjs:
      saveobjH = self.cacheFileHandle[savobj]
      if savobj=='base_dataset_instance':
        savelist = self.rdDset(saveobjH['_savelist'])
        for saveobjname in savelist:
          self.config.base.__dict__[saveobjname] = self.rdDset(saveobjH[saveobjname])
        if '_h5list' in saveobjH.keys():
          h5list = self.rdDset(saveobjH['_h5list'])
          for h5objname in h5list:
            self.config.base.__dict__[h5objname] = self.mkH5list(saveobjH[h5objname])
      else:
        if not hasattr(self.config.base,savobj):
          if savobj in detectors:
            self.config.base.detectors.append(savobj)
            self.config.base._initdet(savobj)
          else:
            self.config.base.__dict__[savobj] = dropData()
        savelist = self.rdDset(saveobjH['_savelist'])
        for saveobjname in savelist:
          self.config.base.__dict__[savobj].__dict__[saveobjname] = self.rdDset(saveobjH[saveobjname])
        if '_h5list' in saveobjH.keys():
          h5list = self.rdDset(saveobjH['_h5list'])
          for h5objname in h5list:
            self.config.base.__dict__[savobj].__dict__[h5objname]  = self.mkH5list(saveobjH[h5objname])

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

    if type(data) is list or ixppyList:
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
    if type(data) in [np.ndarray,int,float,str,
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
      if tdtype == 'list':
        listkeys = [key for key in allkeys if '#' in key]
        listkeys.sort()
        data = []
        for key in listkeys:
          data.append(self.rdDset(rootgroup[key]))
      
      if tdtype == 'dict':
        dictkeys = [key for key in allkeys if '_dtype' not in key]
        data = dict()
        for key in dictkeys:
          data[key] = self.rdDset(rootgroup[key])

    else:
      data = rootgroup.value
      if data is 'empty':
        data = []

    return data

############ CONFIG FILE #############
def _rdConfigurationRaw(fina="ixppy_config"):
  """Reads configuration file that has the description of the different detectors, the default data paths etc."""
  if not os.path.isfile(fina):
    path = os.path.abspath(__file__)
    #print path
    fina = os.path.dirname(path) + '/' + fina    
  filecontent = dict()
  file = open(fina)
  foundlabel = False
  while 1:
    line = file.readline()
    if not line:
      try:
        filecontent[tlabel] = tdat
      except:
        pass
      break

    if foundlabel:
      if line[0] is not '#':
        if line.split():
          tdat.append(line.split())
      else:
        # case when dataset should be finished
        filecontent[tlabel] = tdat
        foundlabel = False

    #label line, keyword to characterize dataset starting next line
    if line[:2]=='#*':
      line = line[2:]
      tlabel = line.split()[0]
      foundlabel = True
      tdat = []
  file.close()
  return filecontent

rdConfigurationRaw = _rdConfigurationRaw

def rdConfiguration(fina="ixppy_config",beamline=None):
  dat = _rdConfigurationRaw(fina)
  home = expanduser("~")
  if os.path.exists(home+'/.ixppyrc'):
    datuser = _rdConfigurationRaw(home+'/.ixppyrc')
    dat = tools.dict_merge(dat,datuser)
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
  cnf["defaultPath"] = _interpreteDefaultPath(dat,cnf['beamline'])
  cnf["scan_step"] = dat['scan_step'][0]
  cnf["cache_directory"] = dat['cache_directory'][0]
  cnf["epics_dset"] = dat['epics_dset'][0][0]
  cnf["epics_cc_dset"] = dat['epics_dset'][0][1]

  return cnf

############ FILE/PATH TOOLS  #############
def _interpreteDefaultPath(confraw,beamline):
  cnf = dict()
  for line in confraw["default_datapath_hosts"]:
    if line[0]==beamline:
      cnf[line[1]] = line[2]
  for line in confraw["default_datapath"]:
    if line[0]==beamline:
      cnf['default'] = line[1]
      break
  return cnf

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
      cnf[alias]['dataset_data'] = dat[detname][n][1]
      cnf[alias]['dataset_time'] = dat[detname][n][2]
      cnf[alias]['dataset_conf'] = dat[detname][n][3]
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

def getProfileLimits(Areadet,step=0,shots=range(10),transpose=False):
  I = Areadet.rdStepData(step,shots)
  tools.nfigure('Select limits')
  pl.imshow(np.squeeze(np.mean(I,axis=0)),interpolation='nearest')
  print 'Select region of interest'
  direction = 'vertical'
  if transpose: direction = 'horizontal'
  lims = np.round(tools.getSpanCoordinates(direction))
  limsdict = dict(projection='vertical range',limits=lims)
  if not hasattr(Areadet,'profileLimits'):
    Areadet._add_saved_datafield('profileLimits',[])
  Areadet.profileLimits.append(limsdict)
  return limsdict

############ TIMING TOOL #############

def TTextractFromRun(run,opalNo=0,referenceCode=162,refthreshold=True,laseroffCode=67,lmonthreshold=True,outputfile=None,filter='standard'):
  d = dataset(run)
  d.config.cache_filename = outputfile
  try:
    d.config.cache.delete()
  except:
    print "no cache file found"
  refdat  = getattr(d.eventCode,'code'+str(referenceCode))
  ldat = getattr(d.eventCode,'code'+str(laseroffCode))     
  TTextract(d.__dict__['opal'+str(opalNo)],
            refmon=refdat,lmon=ldat,       
            refthreshold=refthreshold,
            lmonthreshold=lmonthreshold,filter=filter)
  d.save()

def TTstandardfilter():
  path = os.path.abspath(__file__)
  pathname = os.path.dirname(path)
  weights = np.loadtxt(pathname+'/TTfilt_standard.dat')
  filt = dict()
  filt['weights']=weights
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

def TTapplyFilter(data,filtsettings,plotOutput=False,polysettings=None,erfsettings=None,saveplots=False):
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
      mpr = f.argmax()
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



def TTextractProfiles(Areadet, refmon=None, refthreshold=.1, lmon=None,lmonthreshold=0.05,Nxoff=3, profileLimits=None, profileLimitsRef=None,transpose=False, Nmax=None, steps=None,calibcycles=None,append_to_det=True,saveoffs=False,dataset_name_traces='TTtraces'):
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
  profiles = []
  if not profileLimits:
    profileLimits=Areadet.profileLimits[-1]
  if (lmon is None):
    detChunks = Areadet.chunks(pointdets=dict(xI0=refmon),Nmax=Nmax,steps=steps)
  #elif type(lmon) is tuple:
    #detChunks = Areadet.chunks(pointdets=dict(xI0=refmon),Nmax=Nmax,steps=steps)
  else:
    detChunks = Areadet.chunks(pointdets=dict(xI0=refmon,laser=lmon),Nmax=Nmax,steps=steps)
  sigindices  = []

  allcctraces = []
  if (calibcycles is None):
    calibs = range(len(detChunks))
  else:
    calibs = tools.iterfy(calibcycles)
  Ncalbs = len(clibs)
  for ccNO in calibs:
      print 'starting CalibCycle %d of %d' %(ccNO,Ncalibs)
      cc = detChunks[ccNO]
      tchunksig    = []
      tchunkoff    = []
      tchunktime    = []
      for ch in cc:
          print "Started working on a chunk... (%d shots)" % len(ch.xI0)
          # read all data in chunk
          tdat = ch.data
          ttime = ch.time

          # create average images (fallback off image if not off are found ...)
          profAvAll =  None; # will be calculated later only if needed

          # make profiles of all
          prof    = extractProfilesFromData(tdat,profileLimits)
          if (profileLimitsRef is not None):
              profRef = extractProfilesFromData(tdat,profileLimitsRef)

          # find which shots are reference
          if type(refthreshold) is bool:
            xrayoff = ch.xI0==refthreshold
          else:
            xrayoff = ch.xI0<refthreshold
          if (lmon is None):
            laseroff = np.zeros_like(ch.xI0,dtype=np.bool)
          else:
            if type(lmonthreshold) is bool:
              laseroff = ch.laser==lmonthreshold
            else:
              laseroff = ch.laser<lmonthreshold
          #refind = (xrayoff)&(~laseroff)
          xoffinds = xrayoff.nonzero()[0]
          loffinds = laseroff.nonzero()[0]
          loninds = (~laseroff).nonzero()[0]
          Nshots   = len(ch.xI0)
          # find which shots are signal
          #sigind = (~xrayoff)&(~laseroff)
          
          #correct all signal traces for individual reference
          chsig = []
          choff = []
          chtime = []
          # DIRTY HACK
          for shot in loninds:
          #for shot in range(len(loninds)):
              # find closest Nxoff references and 
              # average as many xoff images as requested (Nxoff)
              if (xoffinds.size == 0):
                  if (profAvAll is None):
                    # create average images (fallback off image if not off are found ...)
                    AvAll     = np.array([np.mean(tdat,axis=0)])
                    profAvAll = extractProfilesFromData(AvAll,profileLimits)[0]
                  poff = profAvAll
              else:
                  temp = (shot<xoffinds).argmax()
                  m=int(temp-float(Nxoff)/2); M = int(temp+float(Nxoff/2))
                  if (m<0): m=0
                  if M>Nshots: M=Nshots
                  # indeces of off images to average
                  xoffinds_toav  = xoffinds[ m:M ]
                  # calc signal
                  poff = np.mean(prof[xoffinds_toav,:],axis=0)
              p    = prof[shot,:]
              time_shot = ttime[shot]

              if (profileLimitsRef is None):
                tsig = (p-poff)/poff
                toff = poff
              else:
                pref    = profRef[shot,:]
                prefoff = np.mean(profRef[xoffinds_toav,:],axis=0)
                tsig = (p-poff*pref/prefoff)/poff
                toff = poff
              chsig.append(tsig)
              choff.append(toff)
              chtime.append(time_shot)

          # make ndarray from all signal traces in chunk and plot it
          if chsig:
            chsig = np.vstack(chsig)
          if choff:
            choff = np.vstack(choff)
          if chtime:
            chtime = np.hstack(chtime)
          tchunksig.append(chsig)
          tchunkoff.append(choff)
          tchunktime.append(chtime)

      if tchunksig:
        tchunksig = np.vstack(tchunksig)
        tchunkoff = np.vstack(tchunkoff)
        tchunktime = np.hstack(tchunktime)
      #allcctraces.append(tchunksig)
      if append_to_det:
        Areadet.append_h5_datalist(dataset_name_traces,tchunksig)
        Areadet.append_h5_datalist(dataset_name_traces+'_time',tchunktime)
        if saveoffs:
          print "off"
          Areadet.append_h5_datalist('TTofftraces',tchunkoff)

      else:
        profiles.append(tchunksig)

      print 'ending CalibCycle %d' %(ccNO)
      #figure(1)
      #imshow(tchunksig,interpolation='nearest')
      #axis('normal')
      #axis('tight')
      #draw()
  if not append_to_det:
    return profiles

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

def extractProfilesFromData(data,profileLimits):
  cameraoffset = 32.
  if profileLimits["projection"]=="vertical range":
    profiles = np.mean(data[:,profileLimits['limits'][0]:profileLimits['limits'][1],:],axis=1)-float(cameraoffset)
  return profiles

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

class dropObject:
  pass

class dropData:
  def add_saved_datafield(self,name,data):
    """adds to a list of datasets to be saved. For faster data recovery and for saving custom analysis progresses (e.g. reduced data from pixel detectors)."""
    if not "_savelist" in self.__dict__:
      self._savelist = []
    self._savelist.append(name)
    self.__dict__[name] = data
  pass

def iterfy(iterable):
    if isinstance(iterable, basestring):
        iterable = [iterable]
    try:
        iter(iterable)
    except TypeError:
        iterable = [iterable]
    return iterable

def isiter(iterable):
    if isinstance(iterable, basestring):
        isv = False
    try:
        iter(iterable)
        isv = True
    except TypeError:
        isv = False
    return isv

def iterdepth(iterable):
  """only for lists/arrays, only along first element"""
  if isiter(iterable):
    N = 0
    iter = True
    while iter:
      N+=1
      iter = eval('isiter(iterable'+(N)*'[0]'+')')
  else: 
    N=0
  return N

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


#def digitizeDataSingle(data,bins,prcentile=.9):
  #"""Binning data"""
  #dims = len(bins)
  #issinglestep = True
  #if len(data) is not dims:
    #if dims==1:
      #data = [data]
    #else:

      #print "binData: number of data structures does not fit number of bin arrays"
      #return

  #abins = []
  #aN = []
  #for dim in range(dims):
    #tbins = bins[dim]
    #if type(tbins)==int or type(tbins)==float:
      ## automatic binning on basis of all cc data
      #adata = data[dim]

      #if tbins<0:
        ## automatic limits
        #if type(tbins)==int:
          ## number of bins, smart finding of borders, prcentile as kwarg
          #pass
        #elif type(tbins)==float:
          ## bin size, automatic finding of borders, prcentile as kwarg 
          #pass
        

      #elif tbins==0:
        ## graphical input of bin edges
        #pass
      #else:
        ## tbins>0
        #if type(tbins) is int:
        ## number of bins given
          #cbins = np.linspace(min(adata),max(adata),tbins+1)

        #elif type(tbins) is float:
        ## bin size given, full range
          #cbins = np.arang(min(adata),max(adata),tbins)
        #pass

    #else:
      ## bin edges given
        #cbins = tbins

    #abins.append(cbins)
    #tdata = data[dim]
    #tN = np.digitize(adata,cbins)
    #aN.append(tN)

  ## initialize index matrix
  #stepbinmatshape = [len(bb)+1 for bb in abins]
  #binmatshape = [len(aN[0])]
  #binmatshape.extend(stepbinmatshape)
  #binmat = np.empty(binmatshape,dtype=np.object_)
  #binmat.fill([])
  #binmat = np.frompyfunc(list,1,1)(binmat)

  #for evNo in range(len(aN[0])):
    #tind = [stepNo]
    #tind.extend([dd[stepNo][evNo] for dd in aN])
    ##print tind
    ##raw_input()
    #binmat[tuple(tind)].append(evNo)

  #if issinglestep:
    #binmat = binmat[0]

  #return binmat,abins

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
        cbins = np.linspace(lims[0],lims[1],np.abs(nbins))

      else:
        # tbins>0
        if type(tbins) is int:
        # number of bins given
          cbins = np.linspace(min(adata),max(adata),tbins+1)

        elif type(tbins) is float:
        # bin size given, full range
          cbins = np.arang(min(adata),max(adata),tbins)
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
  return np.array(bindat)

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
    
  return np.array(dmean),np.array(dstd),np.array(dmedian),np.array(dmad),np.array(dsum),np.array(dN)

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



def calc_polyFitPar(scanvec,detector,monitor,binning=None):
  pass

def calc_weighted_ratio(detector,monitor,binning=None):
  pass

def plotScan(data, monitorIPM='ipm3', detector='diodeU', detectorfieldname='channel0', monitorFiltLims=None, detectorFiltLims=None, scanvec=[], binning=None,binningFiltLims=1,figName=None, oversampling=1,binningCalib=None,centerbinning=False):
  """ NB the binning assumes that the size of elements matches the dataset after timestamp cleanup"""
  #de=bug
  #if isinstance(data,ixppy.dataset):
  if type(data) is not str or not tuple:
    d=data
  else:
    datasets = [monitorIPM,detector]
    if binning:
      if type(binning) is tuple:
        datasets.append(binning[0])
    d = dataset(data,datasets)
    d.filtTimestamps()
  M = d.__dict__[monitorIPM].sum
  I = d.__dict__[detector].__dict__['__'+detectorfieldname]()
  if binning:
    if type(binning) is tuple:
      B = d.__dict__[binning[0]].__dict__['__'+binning[1]]()
    else:
      B = binning 
  if not monitorFiltLims:
    parameterFilt(M,d,name='Amonitor',figName='plotScan monitor filter')
  else:
    parameterFilt(M,d,name='Amonitor',lims=monitorFiltLims)

  if detectorFiltLims:
    if  detectorFiltLims==1:
      parameterFilt(I,d,name='Bdetector',figName='plotScan detector filter')
    elif len(detectorFiltLims)==2:
      parameterFilt(I,d,name='Bdetector',lims=detectorFiltLims)

  if len(scanvec)==0:
    scanvec = d.scanVec
  if binning:
    if binningFiltLims:
      if  binningFiltLims==1:
        binFlt = parameterFilt(B,d,name='Cbinning',figName='plotScan binning filter')
      elif len(binningFiltLims)==2:
        parameterFilt(B,d,name='Cbinning',lims=binningFiltLims)
      B = [ b[~flt] for b,flt in zip(B,d._mergefilters()) ]
    if binningCalib is not None:
      B = [np.polyval(binningCalib,b) for b in B]
    cB = np.mean(np.hstack(B))
    print cB
    if centerbinning:
      Bc = [b-cB for b in B]
    else:
      Bc = B
    binscanvec = [sv+bc for sv,bc in zip(scanvec,Bc)]
    bins = tools.oversample(scanvec,oversampling)
    binsz = np.mean(np.diff(bins))
    edges = np.hstack([bins-binsz/2,bins[-1]+binsz/2])
    indxs = np.digitize(np.hstack(binscanvec),edges)
    #de=bug


    M = d.__dict__[monitorIPM].sum
    I = d.__dict__[detector].__dict__['__'+detectorfieldname]()
    Inorm = np.bincount(indxs,weights=np.hstack(I),minlength=len(edges)+1) / np.bincount(indxs,weights=np.hstack(M),minlength=len(edges)+1)
    Inorm = Inorm[1:-1]
    scanvec = bins
    #de=bug

  else:
    M = d.__dict__[monitorIPM].sum
    I = d.__dict__[detector].__dict__['__'+detectorfieldname]()
    Inorm = []
    for m,i in zip(M,I):
      if not len(m)==0 and not len(i)==0:
        Inorm.append(sum(i)/sum(m))
      else:
        Inorm.append(np.nan)
    Inorm = np.array(Inorm)
  
  if not figName:
    figName = 'plotScan figure'

  tools.nfigure(figName)
  #de=bug
  try:
    pl.plot(scanvec,Inorm,'.-')
  except:
    pl.plot(scanvec[:-1],Inorm,'.-')

  pl.xlabel(d.scanMot)
  pl.ylabel('$\sum$'+detector+'.'+detectorfieldname+' / $\sum$'+monitorIPM+'sum')
  pl.draw()
  return scanvec,Inorm









###### Filtering events ############

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

def parameterCorrFilt(par0,par1,dataset=None,name=None,lims=None,graphicalInput=True,scanstep=None,figName=None,ratio=False):
  print dataset
  par0 = iterfy(par0)
  par1 = iterfy(par1)
  if dataset:
    dataset.filtTimestamps()

  if not scanstep:
    dat0 = np.hstack([spar for spar in par0])
    dat1 = np.hstack([spar for spar in par1])
  else:
    dat0 = par0[scanstep]
    dat1 = par1[scanstep]
  if not figName:
    figName = 'Select filter limits'

  if graphicalInput and not lims:
    tools.nfigure(figName)
    pl.clf()
    if ratio:
      pl.plot(dat0,dat1/dat0,'.k',ms=1)
    else:
      pl.plot(dat0,dat1,'.k',ms=1)
    lims = tools.getRectangleCoordinates()
    lims = list(np.reshape(lims,[2,-1]))



  tfilt0 = []
  tfilt1 = []
  for tpar0,tpar1 in zip(par0,par1):
    if ratio:
      print np.shape(~((tpar0>lims[0][0])&(tpar0<lims[0][1])))
      print np.shape(~((tpar1/tpar0>lims[1][0])&(tpar1/tpar0<lims[1][1])))
      tfilt0.append(~((tpar0>lims[0][0])&(tpar0<lims[0][1])))
      tfilt1.append(~((tpar1/tpar0>lims[1][0])&(tpar1/tpar0<lims[1][1])))
    else:
      print np.shape(~((tpar0>lims[0][0])&(tpar0<lims[0][1])))
      print np.shape(~((tpar1>lims[1][0])&(tpar1<lims[1][1])))
      tfilt0.append(~((tpar0>lims[0][0])&(tpar0<lims[0][1])))
      tfilt1.append(~((tpar1>lims[1][0])&(tpar1<lims[1][1])))

  if dataset:
    #if not name:
      #name = 'unnnamed'
    #filtname = '_filt_'+name
    #dataset.__dict__[filtname] = tfilt
    #dataset._initfilter()
    #return lims
    #dataset.filter = tfilt0
    #print 'tfilt0 done!!!'
    tfilt = [tf1|tf0 for tf1,tf0 in zip(tfilt1,tfilt0)]
    dataset.filter = tfilt
    dataset._initfilter()
    return lims
  else:
    return tfilt0,tfilt1,lims


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
  d = [np.array([c0,c1,c2,c3]).transpose() for c0,c1,c2,c3 in zip(ipmdet._channel0,ipmdet._channel1,ipmdet._channel2,ipmdet._channel3)]
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

def applyOperator(optr,a,b,isreverse=False):
  a = iterfy(a)
  b = iterfy(b)

  res = ixppyList()
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

def _applyFun(func,a):
  res = ixppyList()
  for ta in a:
    res.append(func(ta))
  return res


class ixppyList(list):
  #self.config = config
  #def __new__(self,dset=None):
    #self.dset = dset
    #return self
  def ravel(self):
    return np.hstack(self)
  R = property(ravel)
  #def setDsetFilt(self,**kwargs):
    #parameterFilt(self,

  #self.__dict__['__'+fieldname] = partial(self._compress_name,fieldname)
  #setattr(self.__class__,fieldname,property(self.__dict__['__'+fieldname]))
  def __add__(self,other):
    return applyOperator(operator.add,self,other)
  def __radd__(self,other):
    return applyOperator(operator.add,self,other)
  def __mul__(self,other):
    return applyOperator(operator.mul,self,other)
  def __rmul__(self,other):
    return applyOperator(operator.mul,self,other)
  def __div__(self,other):
    return applyOperator(operator.div,self,other)
  def __rdiv__(self,other):
    return applyOperator(operator.div,self,other,isreverse=True)
  def __truediv__(self,other):
    return applyOperator(operator.truediv,self,other)
  def __rtruediv__(self,other):
    return applyOperator(operator.truediv,self,other,isreverse=True)
  def __floordiv__(self,other):
    return applyOperator(operator.floordiv,self,other)
  def __rfloordiv__(self,other):
    return applyOperator(operator.floordiv,self,other,isreverse=True)
  def __mod__(self,other):
    return applyOperator(operator.mod,self,other)
  def __rmod__(self,other):
    return applyOperator(operator.mod,self,other,isreverse=True)
  def __sub__(self,other):
    return applyOperator(operator.sub,self,other)
  def __rsub__(self,other):
    return applyOperator(operator.sub,self,other,isreverse=True)
  def __pow__(self,other):
    return applyOperator(operator.pow,self,other)
  def __rpow__(self,other):
    return applyOperator(operator.pow,self,other,isreverse=True)

  def __and__(self,other):
    return applyOperator(operator.and_,self,other)
  def __rand__(self,other):
    return applyOperator(operator.and_,self,other)
  def __or__(self,other):
    return applyOperator(operator.or_,self,other)
  def __ror__(self,other):
    return applyOperator(operator.or_,self,other)
  def __xor__(self,other):
    return applyOperator(operator.xor,self,other)
  def __rxor__(self,other):
    return applyOperator(operator.xor,self,other)
  def __le__(self,other):
    return applyOperator(operator.le,self,other)
  def __lt__(self,other):
    return applyOperator(operator.lt,self,other)
  def __eq__(self,other):
    return applyOperator(operator.eq,self,other)
  def __ne__(self,other):
    return applyOperator(operator.ne,self,other)
  def __ge__(self,other):
    return applyOperator(operator.ge,self,other)
  def __gt__(self,other):
    return applyOperator(operator.gt,self,other)


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
      cnf['areaDet'][detname]['dataset_conf'] = 'dummy'
      cnf['areaDet'][detname]['dataset_time'] = dset['dset_time'].split(ccn)[1]
      cnf['areaDet'][detname]['dataset_data'] = dset['dset_data'].split(ccn)[1]
    else:
      cnf['pointDet'][detname] = dict()
      cnf['pointDet'][detname]['dataset_conf'] = 'dummy'
      cnf['pointDet'][detname]['dataset_time'] = dset['dset_time'].split(ccn)[1]
      cnf['pointDet'][detname]['dataset_data'] = dset['dset_data'].split(ccn)[1]
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






############## GENERAL STUFF ###############

beamlines = ['amo','sxr','xpp','xcs','cxi','mec','mob']

cnfFiles = dict()
for bl in beamlines:
  cnfFiles[bl] = rdConfiguration(beamline=bl)
  exec(bl+'=_Lcls_beamline(cnfFiles[bl])')

cnfFile = rdConfiguration()
beamline = cnfFile['beamline']
# initialize pointdet readers
for det in cnfFile['pointDet'].keys():
  exec('rd'+det[0].upper()+det[1:]+'AllData = partial(_rdDetAllData,det=det)')

# initialize areadet readers
for det in cnfFile['areaDet'].keys():
  exec('rd'+det[0].upper()+det[1:]+'StepData = partial(_rdDetStepData,det=det)')

point_detectors = cnfFile['pointDet'].keys()
area_detectors = cnfFile['areaDet'].keys()
detectors = point_detectors
detectors.extend(area_detectors)
