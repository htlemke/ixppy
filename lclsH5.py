from __future__ import print_function
import h5py
import tools
import toolsHdf5 as tH5
import time
import re
import numpy as np
from toolsLog import logbook
#import lclsH5methods as mmm
from toolsHdf5 import datasetRead as h5r
from functools import partial
import os

class lclsH5(object):
  """Class to be defined once to do all hdf5 specific stuff for lcls H5 files."""
  def __init__(self,files,cnf):
    self.fileNames = files
    self.detectors = []
    self.cnf = cnf

  def checkFiles(self):
    self.fileHandles = tools.iterate(self.fileNames,tH5.openOrCreateFile,"r",driver='sec2')

  def close(self):
    for fh in self.fileHandles:
      fh.close()

  def findDetectors(self,detectors=None,exclude=None): 
    # Detectors
    # TODO
    # strategy in 2 paths:
    # (1) detectors are given as alias --> try to read with cached datasets/adresses --> if fails use dataset finder ONLY for defined aliases.
    # (2) no detectors are given --> find all detectors for all aliases and get dataset name if alias not existing.
    #
    if (detectors != "parse"):
      # _findDetectors tries to match datasets found in files with mnemonic given in config
      t0 = time.time()
      logbook("Finding data in hdf5 file ...",end="")
      self.findDatasets(detectors=detectors,exclude=exclude)
      logbook(" ... done (%.1f) ms" % ((time.time()-t0)*1e3),time=False)
    else:
      logbook("Starting to look in the file")
      # parsing look in the file for dataset ...
      # use first (data or cached) file to find detectors to use
      h = self.fileHandles[0]
      try:
        cnf = parseToCnf(h)
        self.cnf = tools.dictMerge(self.cnf,cnf)
        self.areaDet  = cnf['areaDet'].keys()
        self.pointDet = cnf['pointDet'].keys()
        self.detectors = cnf['areaDet'].keys()+cnf['pointDet'].keys()
      except KeyError:
        logbook("Failed to find detectors in ", h.filename)

  def findDatasets(self,detectors=None,exclude=None):
    """finds datasets from a cnf that contains aliases, if no aliases are defined the file is parsed and the hdf5 names are returned as names.
    
    Finds detectors in hdf5 file matching with mnemonic given in config file;
    the matching mnemonic names are as dictionaries (self.pointDet and self.areaDet)
    The 
    """
    subSelection = detectors
    if (subSelection==[]) or (subSelection is None):
      subSelection = self.cnf["pointDet"].keys() + self.cnf["areaDet"].keys()

    if exclude is not None:
      exclude = tools.iterfy(exclude)
      for tex in exclude:
        while True:
          try:
            subSelection.remove(tex)
            continue
          except:
            break
    h = self.fileHandles[0]
   
    # Getting all Detector path strings in CCs and config
    try:
      # try to use only CalibCycle0
      # bad for MEC as some calib cycles don't contain amything... look for longest dataset for now, later look in all

      base = "Configure:0000/Run:0000/"
      bases = h[base].keys()
      lens = np.array([len(h[base][key].keys()) for key in bases])
      base = base + bases[lens.argmax()] +'/'
      h5names = tH5.getDataset_hack(h[base])
      #h5names = [base+x for x in h5names]
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

    

    #raise NotImplementedError('Use the source, luke!')
    ret = {}
    ## *** start EpicsPV *** #
    ## look for epics name
    #epicsFound=False
    #if ("epics_dset" in self.cnf):
      #epicsMne = self.cnf["epics_dset"][0]
      #epicsReg = self.cnf["epics_dset"][1]
      #epicsH5Names=[x for x in h5names if (x.find(epicsReg)>-1)]
      ## common Epics path:
      #ntemp = min([len(x.split("/")) for x in epicsH5Names])
      #epicsCommon = "/".join(epicsH5Names[0].split("/")[0:ntemp])
      ## epics var
      #self._epicsPaths = {}
      #for d in h[epicsCommon]:
        #dpath = d
        #d = d.replace(':','_')
        #d = d.replace('-','_')
        #d = d.replace(' ','_')
        #d = d.replace('.','_')
        #mne = "%s.%s" % (epicsMne.split("/")[0],d)
        #self._epicsPaths[mne]={}
        #self._epicsPaths[mne]["data"] = epicsCommon.replace('CalibCycle:0000','CalibCycle:%04d')+"/"+dpath+"/data"
        #self._epicsPaths[mne]["time"] = epicsCommon.replace('CalibCycle:0000','CalibCycle:%04d')+"/"+dpath+"/time"
        #self._epicsPaths[mne]["conf"] = []
      #self._epicsNames = self._epicsPaths.keys()
    #else:
      #self._epicsNames = []
    ## *** stop EpicsPV *** #
    pointDet = self.cnf["pointDet"]
    for (mnemonic,name) in pointDet.iteritems():
      if (mnemonic.find("nops")>-1) and (mnemonic.find("*")>-1):
        continue
      mnemonic = mnemonic.split('_bak')[0]
      # skip if not in the group we want to read
      if mnemonic not in subSelection:
        continue
      nameData = name["data"].replace("*","\S+")
      detDataset = [x for x in h5names if (re.search(nameData,x) is not None)]
      nameConf = name["conf"].replace("*","\S+")
      try:
        detConf    = [x for x in h5confs if (re.search(nameConf,x) is not None)]
      except:
              detConf=[]
      data = [x for x in detDataset if x[-5:]=="/data" or x[-8:]=="/evrData" or x[-13:]=="/channelValue"]
      time = [x for x in detDataset if x[-5:]=="/time"]
      if ( (len(data) != 0) and (len(time) != 0) ):
        ret[mnemonic] = {}
        #ret[mnemonic]["data"] = data[0].replace('CalibCycle:0000','CalibCycle:%04d')
        #ret[mnemonic]["time"] = time[0].replace('CalibCycle:0000','CalibCycle:%04d')
        ret[mnemonic]["data"] = [replaceCalibCycleString(tdat) for tdat in data]
        ret[mnemonic]["time"] = [replaceCalibCycleString(ttim) for ttim in time]
        if len(detConf)>0:
          ret[mnemonic]["conf"] = detConf[0]
    self._pointDetPaths = ret
    self.pointDetNames = ret.keys()



    areaDet = self.cnf["areaDet"]
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
      #raise NotImplementedError('Use the source, luke!')
      if ( (len(data) != 0) and (len(time) !=0) ):
        ret[mnemonic] = {}
        ret[mnemonic]["data"] = [replaceCalibCycleString(tdat) for tdat in data]
        ret[mnemonic]["time"] = [replaceCalibCycleString(ttim) for ttim in time]
        ret[mnemonic]["conf"] = conf
    self._areaDetPaths = ret
    self.areaDetNames = ret.keys()
    self._detectorsPaths = tools.dictMerge(self._pointDetPaths,self._areaDetPaths)
    self.detectorsNames = self.pointDetNames + self.areaDetNames
    # *** start scan variables *** #
    temp = dict()
    if (len(self.cnf["scan_step"])>0):
      for scan_var in self.cnf["scan_step"]:
        mne,reg = scan_var
        reg  = reg.replace("*","\S+")
        data = [x for x in h5names if (re.search(reg,x) is not None)]

    if not data==[]:
      path = replaceCalibCycleString(data[0])
      obj = scanVar(self.fileHandles,mne,path)
      temp[mne] = obj
    self.scanVars = temp
    # *** stop scan variables *** #
    return

  def initDetectors(self,detectors=None):
    detectors = self.detectorsNames
    # define detectors
    self.detectors = {}

    # START POINT DETECTORS
    t0 = time.time()
    logbook("defining pointDet (with memory cache) ...",end="")
    #if (rdPointDetectorsImmediately):
      #print " (pre-reading all) ",
    for dname in self.pointDetNames:
      #TODO: here dtector dependent modules from a folder are to be used in special cases, like a plugin. Not working because of importing issues. Commented out for now.
      if dname in pluginNames:
        tdclass = eval('detector_'+dname)
      else:
        tdclass = detector
      det = tdclass(self.fileHandles,dname,self._detectorsPaths[dname],useMemoryCache=True,isPointDet=True)
      self.detectors[dname] =det

      #if (rdPointDetectorsImmediately):
        #for i in range(len(self.fileHandles)):
          #det.readData(stepSlice=range(det._numOfScanSteps[i]),fileSlice=i)
      #tools.addToObj( self,dname,det )

    logbook(" ... done (%.1f) ms, %d detectors" % 
      ((time.time()-t0)*1e3,len(self.pointDetNames)),time=False)
    
    # DONE POINT DETECTORS

    # START AREA DETECTORS
    t0 = time.time()
    logbook("defining areaDet (without memory cache) ...",end="")
    for dname in self.areaDetNames:
      det = detector(self.fileHandles,dname,self._detectorsPaths[dname],useMemoryCache=False)
      self.detectors[dname] =det
      #tools.addToObj( self,dname,det )
    logbook(" ... done (%.1f) ms, %d detectors" %
      ((time.time()-t0)*1e3,len(self.areaDetNames)),time=False)
    # DONE AREA DETECTORS
    # consistency check (when daq crashes, last calibs may not have all detectors)
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
    Nc = np.array(Nc)
    NcMin = Nc.min(axis=0)
    NcMax = Nc.max(axis=0)
    for d in tocheck:
      # print out most limiting detector
      if (list(NcMin) == list(d._numOfScanSteps)) and (list(NcMin)!=list(NcMax)):
        logbook("WARNING: Detector/Scan ",d,"is limiting the number of Calybcycle to",str(NcMin),"instead of ",str(NcMax))
      d._numOfScanSteps = list(NcMin)
    self.numOfScanSteps = list(NcMin)
    if len(NcMin) ==1:
      self.numOfScanSteps = self.numOfScanSteps[0]

  
  
class detector(object):
   
  #def __dir__(self):
    #self.__init__heavy()
    #return self.__dict__.keys()

  def __init__(self,fhandle,name,paths,useMemoryCache=True,isauto=False,isPointDet=False):
    self._h5s = tools.iterfy(fhandle)
    self.name = name
    self._paths = paths
    self._useMemoryCache=useMemoryCache
    self._isPointDet = isPointDet
    self._numOfScanStepsFile = self._checkNcalib()
    #if not isauto:
      #self.__init__heavy()
    if self._isPointDet:
      self._initPointDet()
    else:
      self._initAreaDet()
  
  def _getTotalNumOfSteps(self):
    return np.sum(self._numOfScanStepsFile)
  Nsteps = property(_getTotalNumOfSteps)

  def _getFileStep(self,step):
    stepout = step
    fileind = 0
    while stepout >= self._numOfScanStepsFile[fileind]:
      stepout -= self._numOfScanStepsFile[fileind]
      fileind += 1
    return fileind, stepout
  
  def _initPointDet(self):
    if self._isPointDet and self._useMemoryCache:
      if len(self._paths['data'])>1:
        if len(self._paths['data'])>50:
          # this seems the epics case
          logbook("crazy amount of point counters in %s, force full reading initialization with..." %(self.name))
        else:
        # this might be the tt case
          for datapath,timepath in zip(self._paths['data'],self._paths['time']):
            datapath = getPath(datapath)
            timepath = getPath(timepath)
            name = getDetNamefromPath(datapath)
            dat,fields = self._readPointDataGeneral(datapath)
            times = self._readTime(timepath)
            if not fields==[] and not hasattr(self,'fields'):
              self.fields = dict()
            for field,tdat in zip(fields,dat):
              tname = name+'_'+field
              data = [tdat,times]
              self.fields[tname] = data 
        	#memdata(tname,data)
      else:
         # this might be the ipm case
        datapath = self._paths['data'][0]
        timepath = self._paths['time'][0]
        datapath = getPath(datapath)
        timepath = getPath(timepath)
        dat,fields = self._readPointDataGeneral(datapath)
        times = self._readTime(timepath)

        if not fields==[]:
          self.fields = dict()
          for field,tdat in zip(fields,dat):
            tname = field
            data = [tdat,times]
            self.fields[tname] = data 

  def _initAreaDet(self):
    self.time = self._readTime(getPath(self._paths['time'][0]))
    self.readData = self._readDataGeneral
    self._numOfScanStepsFile = self._checkNcalib()
    self.readData = self._readDataGeneral


      
      #######################################################################
      #######################################################################
      #######################################################################
  def _readPointDataGeneral(self,path,field=None,stepSlice=None,shotSlice=None):
    if stepSlice is None:
      stepSlice = range(self.Nsteps)
    stepSlice = tools.iterfy(stepSlice)

    #self.readTime(stepSlice,shotSlice,fileSlice)
    #if not self._existsInSelf("time"):
      #timeStampObj = memdata(self,"timestamp",fileSlice[0])
      ##self._addToSelf("time",timeStampObj)

    outS = []
    for stepNum in stepSlice:
      fileNum,FstepNum = self._getFileStep(stepNum)
      cpath = path % FstepNum # really beautiful, same for python 3 ?
      cpath = getPath(cpath)
      data = h5r(self._h5s[fileNum],cpath)
      try:
        if (shotSlice is None):
          data = data[...]
        else:
                data = data[shotSlice]
      except:
        data = np.array([])

      outS.append(data)

    # find structure with something in
    outSind = 0
    for toutS in outS:
      if len(toutS)>0:
        break
      outSind+=1

    if outS[outSind].dtype.names:
      if not field is None:
        index = ''
        while not field in outS[outSind].dtype.names:
          index = field[-1] + index
          field = field[:-1]
        index = int(index)
        fields = [field]
      else:
        fields = outS[outSind].dtype.names
        index = None
      
      pret = [[dd[tfield] if len(dd)>0 else np.array([]) for dd in outS ] for tfield in fields]
      ret = []
      retfields = []
      for tret,tfield in zip(pret,fields):
        if tret[0].ndim==2:
          noofvecs = np.shape(outS[0][tfield])[1]
          if not index is None:
            indices = [index]
          else:
            indices = range(noofvecs)
          for tindex in indices:
            strfmt = '%0' + '%dd' %(1+np.floor(np.log10(noofvecs)))
            tname = tfield + strfmt %(tindex)
            ret.append([sd[:,tindex] if np.ndim(sd)==2 else np.asarray([]) for sd in tret ])
            retfields.append(tname)
        else:
          ret.append(tret)
          retfields.append(tfield)
    
    return ret,retfields
      
      #######################################################################
      #######################################################################
      #######################################################################

  

  def __init__heavy(self):
      # based on det type point readData to a particular method
      # container for different subfields or run
      self._guessDetType()
      self._numOfScanStepsFile = self._checkNcalib()
      if (self._useMemoryCache):
        self.readData(stepSlice=0,shotSlice=None,fileSlice=0)

  #def __repr__(self):
    #if (self._useMemoryCache):
      #return "`detector` object: %s (w memory cache)" % self.name
    #else:
      #return "`detector` object: %s (w/o memory cache)" % self.name

  # data object (self)manipulation
  def _getFromSelf(self,what):
    """ get data from object for example: d.ipm2.get("_file0.step3.channel") """
    return tools.getFromObj(self,what)
  def _existsInSelf(self,what):
    return tools.existsInObj(self,what)
  def _addToSelf(self,what,value,overWrite=True):
    return tools.addToObj(self,what,value,overWrite=overWrite)

  def _checkNcalib(self):
    numOfScanStepsFile = []
    path = tools.commonPathPrefix(self._paths["data"])
    path = getPath(path)
    for h in self._h5s:
      n=0
      while( tH5.datasetExists(h,path % n) ):
        n+=1
      numOfScanStepsFile.append(n)
    nccs = []
    for h in self._h5s:
      n=0
      while( tH5.datasetExists(h,'/Configure:0000/Run:0000/CalibCycle:%04d/' % n) ):
        n+=1
      nccs.append(n)
    self._numOfScanStepsFile = []
    for n in range(len(self._h5s)):
      if nccs[n] > numOfScanStepsFile[n]:
        logbook("More calibcycle structures than detectors, will lead to empty detector steps...")
        self._numOfScanStepsFile.append(nccs[n])
      else:
        self._numOfScanStepsFile.append(numOfScanStepsFile[n])

    
    return self._numOfScanStepsFile

  def __getitem__(self,x):
    n = self._numOfScanSteps[0]
    if isinstance(x,slice):
      return [self[ii] for ii in xrange(*x.indices(n))]
    else:
      x = tools.iterfy(x)
      if (max(x)>=self._numOfScanSteps[0]):
        raise IndexError
      if (not self._useMemoryCache):
        return self.readData(stepSlice=x)
      if self._existsInSelf("_data"):
        return self._getFromSelf("_data")[x]
      elif self._existsInSelf("value"):
        return self._getFromSelf("value")[x]
      else:
        return tools.getFromObj(self,x)
  
  def _guessDetType(self):
    data_type = self._h5s[0][self._paths["data"] % 0].dtype
    if (data_type == np.dtype([('channel', '<f4', (4,)), ('sum', '<f4'), ('xpos', '<f4'), ('ypos', '<f4')]) ):
      self._type="IPM"
      self.readData = self._readDataIPM
    elif self.name == "eventCode":
      self._setupEventCodes()
      self.readData = self._readDataEVR
    else:
      self._type="General"
      self.readData = self._readDataGeneral

  def _readTime(self,path,stepSlice=None,shotSlice=None):
    if stepSlice is None:
      stepSlice = range(self.Nsteps)
    stepSlice = tools.iterfy(stepSlice)
    times = []
    for stepNum in stepSlice:
      fileNum,FstepNum = self._getFileStep(stepNum)
      tpath = path % FstepNum
      time = h5r(self._h5s[fileNum],tpath)
      try:  
        if (shotSlice is None):
          time = time[...]
        else:
          time = time[shotSlice]
      except:
        time = np.array([])

      #raise NotImplementedError('Use the source, luke!')
      times.append(time)
    return times

    #if (shotSlice is None):
      #return times
    #else:
      #return time[shotSlice]

  def _setupEventCodes(self):
    h=self._h5s[0]
    codes = h[self._paths["conf"]][...]["code"]
    self._addToSelf("codes",codes)
    #if (EvrFound):

  def _readDataEVR(self,stepSlice=0,shotSlice=None,fileSlice=0):
    fileSlice = tools.iterfy(fileSlice)
    stepSlice = tools.iterfy(stepSlice)
    for fileNum in fileSlice:
      for stepNum in stepSlice:
        data= self._readDataGeneral(stepSlice=stepNum,shotSlice=shotSlice,
          fileSlice=fileSlice)
        timeStampObj = self._getFromSelf("time")
        addr = address(fileNum,stepNum,"shotSlice")
        addrCode = address(fileNum,stepNum,"code%d"%self.codes[0])
        # if not read or we are asking shots outside the range of read values...
        if ( (not self._existsInSelf(addrCode)) or (len(tools.iterDiff( self._getFromSelf(addr), shotSlice) )==0) ):
          for code in self.codes:
            nshots = len(data)
            temp = np.zeros(nshots,dtype=np.bool)
            for nshot in range(nshots):
              if code in data[nshot][0]["eventCode"]:
                temp[nshot] = True
            addr = address(fileNum,stepNum,"code%d"%code)
            self._addToSelf(addr,temp)
            dataObj = memdata(self,"code%d"%code,fileNum,timeStampObj)
            self._addToSelf("code%d"%code,dataObj)

  def _readDataIPM(self,stepSlice=0,shotSlice=None,fileSlice=0,fieldName=None):
    # further splits channels
    fileSlice = tools.iterfy(fileSlice)
    stepSlice = tools.iterfy(stepSlice)
    for fileNum in fileSlice:
      data= self._readDataGeneral(stepSlice=stepSlice,shotSlice=shotSlice,
        fileSlice=fileNum,fieldName=fieldName)
      timeStampObj = self._getFromSelf("time")
      for stepNum in stepSlice:
        addr = address(fileNum,stepNum,"shotSlice")
        addrChannel = address(fileNum,stepNum,"channel0")
        # if not read or we are asking shots outside the range of read values...
        if ( (not self._existsInSelf(addrChannel)) or (len(tools.iterDiff( self._getFromSelf(addr), shotSlice) )==0) ):
          # if only one step is read it does not return list ...
          if isinstance(data,list):
            channel_data = data[stepNum]["channel"]
          else:
            channel_data = data["channel"]
          for i in range(4):
            addr = address(fileNum,stepNum,"channel%d"%i)
            self._addToSelf(addr,channel_data[:,i])
            dataObj = memdata(self,"channel%d"%i,fileNum,timeStampObj)
            self._addToSelf("channel%d"%i,dataObj)

  def _readDataGeneral(self,stepSlice=0,shotSlice=None,fileSlice=0,fieldName=None):
    fileSlice = tools.iterfy(fileSlice)
    stepSlice = tools.iterfy(stepSlice)
    #self.readTime(stepSlice,shotSlice,fileSlice)
    #if not self._existsInSelf("time"):
      #timeStampObj = memdata(self,"timestamp",fileSlice[0])
      #self._addToSelf("time",timeStampObj)
    outS = []
    # NB: Here is an issue, doesn't make sense like it is right now...
    for stepNum in stepSlice:
      # check if memory chached exists and the already read values contains what we need..
      #print "TODO: add slice1 in slice2 ..., the current one does not work with None ..."
      #print "r1=range(5,100); r2=range(10,20); [idx for (idx,x) in zip(range(len(r1)),r1) if x in r2]"
      fileNum,FstepNum = self._getFileStep(stepNum)
      addr = address(fileNum,stepNum,"shotSlice")
      if (not self._existsInSelf(addr)) or (len(tools.iterDiff( self._getFromSelf(addr), shotSlice) )==0):

        path = self._paths["data"][0] % FstepNum
        path = getPath(path)
        data = h5r(self._h5s[fileNum],path)
        if (shotSlice is None):
          data = data[...]
        else:
          if isinstance(shotSlice,np.ndarray) and shotSlice.dtype is np.dtype(int):
            tshotSlice = np.zeros([data.len()],dtype=bool)
            tshotSlice[shotSlice]=True
            shotSlice=tshotSlice

          data = data[shotSlice]
      else:
          data = self._getFromSelf(address(fileNum,stepNum,"_data"))
    # store is asked to use memory cache
      if (self._useMemoryCache):
        # save in .fileNum.stepNum._data
        self._addToSelf(address(fileNum,stepNum,"_data"),data)
        self._addToSelf(address(fileNum,stepNum,"shotSlice"),shotSlice)
      if (isinstance(data.dtype.names,tuple)):
        for fieldname in data.dtype.names:
          self._addToSelf(address(fileNum,stepNum,fieldname),data[fieldname])
          if ( not (fieldname in self.__dict__) ):
            timeStampObj = memdata(self,"timestamp",fileNum)
            dataObj = memdata(self,fieldname,fileNum,timeStampObj)
            self._addToSelf(fieldname,dataObj)
      #else: 
        #timeStampObj = memdata(self,"timestamp",fileNum)
        #dataObj = memdata(self,"_data",fileNum,timeStampObj)
        #tools.addToObj(self,"_data",dataObj)
      outS.append(data)
    return outS

# hack for linked calib cycles
def getPath(name):
  if name[:12] == '/CalibCycle:':
    return '/Configure:0000/Run:0000'+name
  else:
    return name

class scanVar(object):
  def __init__(self,fhandle,name,paths):
    self._h5s = fhandle
    self._name = name
    self._paths = getPath(paths)
    self._checkNcalib()
    self._read()

  def _read(self):
    names = self._h5s[0][getPath(self._paths % 0)]
    if names.shape==():
      return
    names = names['name']
    self.name  = list(names)[0]
    self.names = list(names)
    v = []
    for i in range(len(self._h5s)):
      h=self._h5s[i]
      for c in range(self._numOfScanSteps[i]):
        val = h[self._paths % c]["value"]
        v.append(val)
    self.data = v
    #BUGME
    #for n in names:
      #dhelp = memdata(self,n,0)
      #self.__dict__[n,dhelp]
      #tools.addToObj(self,n,dhelp)

  def __repr__(self):
    return "`scanVar` object: %s" % self._name

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

def parseToCnf(fileHandle):
  cnf = dict()
  cnf['areaDet'] = dict()
  cnf['pointDet'] = dict()
  dsets = findDatasetsInHdf5(fileHandle)

  for dset in dsets:
    detname = tools.varName(os.path.split(os.path.split(dset['dset_data'])[0])[1])
    base = "/Configure:0000/Run:0000/"
    bases = h[base].keys()
    lens = np.array([len(h[base][key].keys()) for key in bases])
    base = base + bases[lens.argmax()] +'/'
    #ccn = '/Configure:0000/Run:0000/CalibCycle:0000/'
    ccn = base

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

def getDetNamefromPath(path,lowerit=False):
  name = 'data'
  while name in ['evrData','data','time','image','array','channelValue']:
    path,name = os.path.split(path)
  name = name.lower()
  name = name.replace(':','_')
  name = name.replace('.','_')
  name = name.replace(' ','_')
  name = name.replace('-','_')
  if lowerit:
    while name.find('_')>-1:
      s,e = name.split('_')
      name = s+e[0].upper()+e[1:]
  return name

def address(fileNum,stepNum,what):
  return "_file%d.step%d.%s" % (fileNum,stepNum,what)

def findDatasetsInHdf5(filehandle):
  rundset = filehandle['Configure:0000']['Run:0000']
  ccname = rundset.keys()[0]
  dsets = crawlforDatasets(rundset[ccname])
  return dsets

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
      elif u'data' in itemnames or u'channelValue' in itemnames:
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
### Henrik: This part I think could go to some plugin directory. I put it here anyway, mainly as I failed to get it to work. The problem I ran into is cross imports, import detector class from plugin module and importing detector wrapper from here didn't work.
pluginNames = ['epics','eventCode']
#pluginNames = []

class detector_epics(detector):
  def _initPointDet(self): 

    for datapath,timepath in zip(self._paths['data'],self._paths['time']):
      name = getDetNamefromPath(datapath)
      

      rdAllFunc = partial(self._readEpicsAllData,datapath,timepath)
      if not hasattr(self,'fields'):
        self.fields = dict()
      self.fields[name] = rdAllFunc

  def _readEpicsAllData(self,datapath,timepath):
    
    stepSlice = range(self.Nsteps)
    stepSlice = tools.iterfy(stepSlice)

    outD = []
    outT = []
    for stepNum in stepSlice:
      fileNum,FstepNum = self._getFileStep(stepNum)
      cpath = datapath % FstepNum # really beautiful, same for python 3 ?
      data = h5r(self._h5s[fileNum],cpath)
      outD.append(data['value'])
      time = data['stamp']
      time.dtype.names = ('seconds','nanoseconds')
      outT.append(time)
    return [outD,outT]

class detector_eventCode(detector):
  def _initPointDet(self):
    detector._initPointDet(self)
    h=self._h5s[0]
    codes = h[self._paths["conf"]][...]["code"]
    self._addToSelf("codes",codes)
    self._getCodeFields()
    #if (EvrFound):

  def _getCodeFields(self):
    dat = self.fields['fifoEvents'][0]
    tim = self.fields['fifoEvents'][1]
    dat = [[evt['eventCode'] for evt in step] for step in dat]
    for code in self.codes:
      name = 'code_'+str(code)
      cd = [np.array([True if code in evt else False for evt in step]) for step in dat]
      self.fields[name] = [cd,tim] 

class detector_evts(detector):
  def _initPointDet(self):
    detector._initPointDet(self)
    h=self._h5s[0]
    codes = h[self._paths["conf"]][...]["code"]
    self._addToSelf("codes",codes)
    self._getCodeFields()
    #if (EvrFound):

  def _getCodeFields(self):
    dat = self.fields['fifoEvents'][0]
    tim = self.fields['fifoEvents'][1]
    dat = [[evt['eventCode'] for evt in step] for step in dat]
    for code in self.codes:
      name = 'code_'+str(code)
      cd = [np.array([True if code in evt else False for evt in step]) for step in dat]
      self.fields[name] = [cd,tim] 

def replaceCalibCycleString(input,digits=4):
  fstr = 'CalibCycle:'
  starti = input.find(fstr)
  replstr = input[starti:starti+len(fstr)+digits]
  input = input.replace(replstr,'CalibCycle:%04d')
  return input


