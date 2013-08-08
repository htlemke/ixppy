""" HDF5 UTILS """
import numpy as np
from toolsLog import logbook
import toolsVarious
import h5py
import os

_datasetsInHdf5File={}

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
  """ read name from hdf5, if it does not exist return None """
  #print "H5 read",name
  if (datasetExists(h5handle,name)):
    return h5handle[name]
  else:
    return None

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

def datasetToNumpy(data,slice=None):
  if slice is None:
    return data[...]
  else:
    return data[slice]
