""" HDF5 UTILS """
import numpy as np
from toolsLog import logbook
import toolsVarious
import h5py
import os
from multiprocessing import Pool


_datasetsInHdf5File={}

class ixppyHDF5(h5py.File):
  def __init__(self,name,mode=None, driver=None, libver=None, userblock_size=None, **kwds):
    self.h = h5py.File.__init__(self,name,mode=mode,driver=driver,libver=libver,\
            userblock_size=userblock_size,**kwds)
  def getNames(self):
    if not hasattr(self,"_names"): self._names = getNames(self)
    return self._names

  def datasetGet(self,name):
    if (name in self):
      return self[name]
    else:
      return None

  def getObj(self):
    if not hasattr(self,"obj"): self.obj = Hdf5ToObj(self)
    return self.obj

  def datasetRead(self,name,chunksize=1):
    dataset = self.dataset(name)
    if dataset is None: return None
    size = dataset.size
    n = dataset.shape[0]
    isMultiProcessUseful = size> (1024*1024)
    if (nChunkSize > 1) and isMultiProcessUseful:
      idx = range(n)
      def f(x):
        return dataset[x]
      p = Pool(); # 16-43 ms overhead
      res = p.map_async(f,idx,chunksize=chunksize)

      


def openOrCreateFile(fname,mode="r"):
  if (os.path.isfile(fname)):
    if (mode == "r"):
      if not os.access(fname,os.R_OK):
        raise IOError("Asked to read %s but it is not possible, check permissions" % fname)
        return None
    elif (mode=="r+") or (mode=="a") or (mode=="w"):
      if not os.access(fname,os.W_OK):
        raise IOError("Asked to read/write %s but it is not possible, check permissions" % fname)
        return None
    h5handle=h5py.File(fname,mode)
    logbook("File %s exists already, opening in %s mode" % (fname,mode))
  else:
    logbook("File %s does not exists, creating it" % (fname))
    h5handle=h5py.File(fname,"w")
  return h5handle

def getNames(h5handle):
  out=[]
  h5handle.visit(out.append)
  return out

def datasetExists(h5handle,name):
  return (name in h5handle)

def datasetWrite(h5handle,name,data):
  """ (Over)write a dataset with data, no shape check is done """
  if datasetExists(h5handle,name): del h5handle[name]
  h5handle[name]=data

def datasetRead(h5handle,name):
  """ read name from hdf5, returning the dataset instance, if it does not exist return None """
  #print "H5 read",name
  if (datasetExists(h5handle,name)):
    return h5handle[name]
  else:
    return None

def getData(h5handle,name,slice=None,chunksize=1):
  dataset = getDataset(h5handle,name)
  return datasetToNumpy(dataset,slice=slice,chunksize=chunksize)

def getDataset(h5handle,reRead=False):
  # caching
  if ((h5handle in _datasetsInHdf5File) and (not reRead)):
    return _datasetsInHdf5File[h5handle]
  out = []
  def func(name,obj):
    if (isinstance,h5py.Dataset):
      out.append(name)
  h5handle.visititems(func)
  globals()["_datasetsInHdf5File"][h5handle]=out
  return out

def getDataset_hack(h5handle,reRead=False):
  # caching
  if ((h5handle in _datasetsInHdf5File) and (not reRead)):
    return _datasetsInHdf5File[h5handle]
  out = []
  def checkkeys(handle):
    for key in handle.keys():
      okey = handle[key]
      if isinstance(okey,h5py.Dataset):
	out.append(okey.name)
      elif isinstance(okey,h5py.Group):
	checkkeys(okey)
  checkkeys(h5handle)

  globals()["_datasetsInHdf5File"][h5handle]=out
  return out

def Hdf5ToObj(h5handle):
  ishdf5group = (isinstance(h5handle,h5py.File)) or (isinstance(h5handle,h5py.Group))
  if not ishdf5group:
    h5handle = h5py.File(h5handle,'r')
  ret = toolsVarious.dropObject()
  for h in h5handle:
    name = h.replace(":","_")
    name = name.replace(".","_")
    if not isinstance(h5handle[h],h5py.Dataset):
      ret._add(name,Hdf5ToObj(h5handle[h]))
    else:
      ret._add(name,h5handle[h])
  return ret

def f(arg):
  dataset,x=arg
  return dataset[x]

def f1(dataset,x):
  return dataset[x]


def datasetToNumpy(dataset,sliceSel=None,chunksize=1):
  size = dataset.size
  n = dataset.shape[0]
  if (sliceSel is None): sliceSel = slice(0,n,1)
  isMultiProcessUseful = size> (1024*1024)
  if (chunksize > 1) and isMultiProcessUseful:
    # subdivide indices in chunksize
    start,stop,step = sliceSel.indices(n)
    nC = int(float(stop-start)/step/chunksize+0.5)
    print nC
    args = []
    for i in range(nC):
      s1 = start+i*(chunksize*step)
      s2 = start+(i+1)*(chunksize*step)
      print i,s1,s2
      args.append( (dataset,slice(s1,s2,step) ) )
    print args
    raw_input("Not working yet, use chunksize = 1")
    p = Pool(); # 16-43 ms overhead
    res = p.map_async(f,args,chunksize=1)
    p.close()
    p.join()
    data = np.asarray(res.get())
  else:
    data = dataset[sliceSel]
  return data
