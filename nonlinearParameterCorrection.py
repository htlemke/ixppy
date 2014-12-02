from ixppy import tools
import pylab as plt
import numpy as np
from scipy import linalg,io
import copy



def corrNonlinGetpars(expar,Imat,order=3,exparWP=0,Iwp=None):
  if Iwp is not None:
    Imat -= Iwp
  pol = np.vander(expar-exparWP,order+1)
  
  if Iwp is not None:
    np.delete(pol,-1,1)
  scale = np.sqrt((pol*pol).sum(axis=0))
  pol /= scale
  comps,resid,rnk,singv = linalg.lstsq(pol,Imat)
  comps = (comps.T/scale).T
  if Iwp is None:
    Iwp = np.array(np.matrix(pol)*np.matrix(comps))
  return comps,exparWP,Iwp

#def corrNonlin(expar,Iuncorr,comps,exparWP,Iwp,):
  #pol = np.vander(expar-exparWP,order+1)
 # 
  #if Iwp is not None:
    #np.delete(pol,-1,1)
  


def getCorrectionFunc(order=5,i0=None,Imat=None,i0_wp=1e6,fraclims_dc=[.9,1.1]):
  """ 
  Getting nonlinear correction factors form a calibration dataset consiting of:
    i0     	array of intensity/parameter values the calibration has been made for
    Imat   	2D array of the corresponding reference patterns, in each row 
  		there is one ravelled array of each intensity bin in i0.
    i0_wp	a working point around which a correction polynomial will be
  		developed for each pixel.
    order	the polynomial order up to which will be deveoped.
    fraclims_dc	relative factor for the i0,Imat data limits which are used to
  		determine the working point location.
  
  Returns corrFunc(i,D), a function that takes a flat array of intensity/
  		parameter values as well as a Matrix D of flattened patterns 
		the correction is to be applied on (rows in D are again corresponding
		to each intensity in i).
  """

  msk = tools.filtvec(i0,i0_wp*np.asarray(fraclims_dc))
  p0 = tools.polyFit(i0[msk],Imat[msk,:],2)
  dc = tools.polyVal(p0,i0_wp)
  comps = tools.polyFit(i0-i0_wp,Imat-dc,order,removeOrders=[0])
  compsder = tools.polyDer(comps)
  c = lambda(i): tools.polyVal(comps,i-np.asarray(tools.iterfy(i0_wp)))+dc
  c_prime = lambda(i): tools.polyVal(compsder,i-np.asarray(tools.iterfy(i0_wp)))
  t = lambda(i): (c_prime(i0_wp).T * (i-i0_wp)).T + dc
  
  cprimeic = c_prime(i0_wp)
  dcorr_const = -cprimeic*i0_wp + c(i0_wp) - t(0) 
  def corrFunc(i,D):
    return (i*cprimeic.T + dcorr_const.T + ((D-c(i))*cprimeic/c_prime(i)).T).T
  return corrFunc



