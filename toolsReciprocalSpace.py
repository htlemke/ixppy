import numpy as np
import toolsVarious
import toolsVecAndMat





def gamdel2Qfibold(gamma,delta,alpha,lam):
  gamma = np.array(toolsVarious.iterfy(gamma))
  delta = np.array(toolsVarious.iterfy(delta))

  shpgam = np.shape(gamma)
  shpdel = np.shape(delta)
  if not shpgam==shpdel:
    print "gamma and delta array must have same shape!"
    return
  gamma = gamma.ravel()
  delta = delta.ravel()
  Qs =  2*np.pi/lam * -np.array(toolsVecAndMat.rotmat3D([0,1,0],alpha)*np.mat([
    np.cos(delta)*np.cos(gamma)-1,
    -np.cos(delta)*np.sin(gamma),
    -np.sin(delta)]))
  Qip = np.sign(Qs[1,:])*np.sqrt(Qs[0,:]**2+Qs[1,:]**2);
  Qop = Qs[2,:]
  Qip = Qip.reshape(shpgam)
  Qop = Qop.reshape(shpgam)
  return Qip,Qop

def gamdel2Qfib(gamma,delta,alpha,lam):
  gamma = np.array(toolsVarious.iterfy(gamma))
  delta = np.array(toolsVarious.iterfy(delta))

  shpgam = np.shape(gamma)
  shpdel = np.shape(delta)
  if not shpgam==shpdel:
    print "gamma and delta array must have same shape!"
    return
  gamma = gamma.ravel()
  delta = delta.ravel()
  Qs =  2*np.pi/lam * np.array((-toolsVecAndMat.rotmat3D([0,1,0],-alpha))*np.mat([
    np.cos(delta)*np.cos(gamma)-1,
    -np.cos(delta)*np.sin(gamma),
    -np.sin(delta)]))
  Qip = np.sign(Qs[1,:])*np.sqrt(Qs[0,:]**2+Qs[1,:]**2);
  Qop = Qs[2,:]
  Qip = Qip.reshape(shpgam)
  Qop = Qop.reshape(shpgam)
  return Qip,Qop

def revertElAz(el,az):
  elO = np.arcsin(np.cos(el)*np.sin(az))                                                               
  azO=-np.arctan(np.sin(el)/(np.cos(el)*np.cos(az)))
  return elO,azO


