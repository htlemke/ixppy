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
#function [tth,chi] = gamdel2tthchi(gam,del)
#%[tth,chi] = gamdel2tthchi(gam,del)
#%H. T. Lemke  2009

#chi = atan(tan(del)./sin(gam));
#tth = acos(cos(gam).*cos(del));


#function [gam,del] = tthechi2gamdel(tthe,chi)

#%Henrik T. Lemke April 2007



#gamdelsym = solve(['tan(del)/sin(gam) = tan(' num2str(chi) ')'],...
    #['cos(' num2str(tthe) ')=cos(gam)*cos(del)'],'gam','del');

#gam = double(gamdelsym.gam(2))
#del = double(gamdelsym.del(2))



#function [gam, del, tth] = Qfib2angles(Qinp,Qoop,lam,alp)
#% all anlges input and output in radians
#% Henrik T. Lemke Oct 2008

#% 
#Qinp = Qinp(:);
#Qoop = Qoop(:);
#% lam  = 2*pi;
#% alpha = 0*pi/180;

#phiE = Qfiber2QphiE(Qinp, Qoop, lam, alp);
#phiE = phiE(:);
#del = asin((Qoop*cos(alp) - Qinp.*sin(phiE)*sin(alp))./(2*pi/lam));
#gam = asin(Qinp.*cos(phiE)./(cos(del)*2*pi/lam));
#tth = acos(cos(gam).*cos(del));


