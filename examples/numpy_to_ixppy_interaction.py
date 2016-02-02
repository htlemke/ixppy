""" Examples showing how to wrap numpy arrays as ixppy objects """
import ixppy
import numpy as np

def example_memdata(d,detName="imp3.sum"):
  det = d.get(detName)
  # create dummy data (of the right shape)
  dummy = [np.random.random(nShots) for nShots in det.lens]
  time  = det.time
  # create memdata
  myWrappedNunmpy = ixppy.memdata(name="myname",input=(dummy,time),scan=d.scan)
  test = myWrappedNunmpy.digitize(np.arange(0,1,0.1))
  return myWrappedNunmpy,test

def example_data(d,detName="ipm3.sum",nQ=100):
  det = d.get(detName)
  # create dummy data (of the right shape); they could be the 1D azimuthal
  # averaging for example
  dummy = [np.random.random( (nShots,nQ) ) for nShots in det.lens]
  time  = det.time
  # wrap the numpy array in a class with a __call__ method
  dummywrap = ixppy.datawrap(dummy)
  myWrappedNunmpy = ixppy.data(name="myname",input=dummywrap,time=time,scan=d.scan)
  return myWrappedNunmpy
 
