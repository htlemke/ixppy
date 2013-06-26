""" HDF5 UTILS """
import numpy as np
from toolsLog import logbook
import h5py
import os

def openOrCreateFile(fname,mode="r+"):
	if (os.path.isfile(fname)):
		h5handle=h5py.File(fname,mode)
		logbook("File %s exists already, opening in %s mode" % (fname,mode))
	else:
		logbook("File %s does not exists, creating it" % (fname))
		h5handle=h5py.File(fname,"w")
	return h5handle

def printTree(h5handle):
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
	if (datasetExists(h5handle,name)):
		return h5handle[name]
	else:
		return None
