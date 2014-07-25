import numpy as np
from toolsConstsAndConv import *
from toolsVarious import *

def getLODCMdelay(E,oldE,ds=.6,ID='Si',hkl=(1,1,1)):
    theold = BraggAngle(ID,hkl,E=oldE)
    thenew = BraggAngle(ID,hkl,E=E)
    pathold = ds/sind(2*theold)-ds/tand(2*theold)
    pathnew = ds/sind(2*thenew)-ds/tand(2*thenew)
    delay = (pathnew-pathold)/c_light()
    return delay


gR = 3.175
gD = 231.303
#    gTheta0 = 15.08219; # original calibration
gTheta0 = 15.18219; # calibration for split and delay
gTheta0 = 14.983; # quick calibration Jul22 4am
gTheta0 = 14.9786; # better calibration Jul27 using the inflection point at data resolution and edges of Ni,Ti and Zr
gTheta0 = 14.9694; # better calibration Jul27 using the inflection point of a splinefit only using Ti edge, the Ti edge was best resolved and is closest to the energy for L333.
gTheta0 = 14.9792; # better calibration Jul27 using the inflection point at (probably) improved resolution using a splinefit.
gSi111dspacing = 3.13556044

def thetaToAlio(theta):
  """ Function that converts theta angle (deg) to alio position (mm) """
  t_rad = (theta-gTheta0)*np.pi/180.
  x = gR*(1/np.cos(t_rad)-1)+gD*np.tan(t_rad)
  return x

def alioToTheta(alio):
  """ Function that converts alio position (mm) to theta angle (deg) """
  return gTheta0 + 180/np.pi*2*np.arctan( (np.sqrt(alio**2+gD**2+2*gR*alio)-gD) / (2*gR+alio) )


def thetaToWavelength(theta):
  """ Function that converts theta angle (deg) to wavelength (A) """
  return 2*gSi111dspacing*np.sin(theta/180*np.pi)

def wavelengthToTheta(wavelength):
  """ Function that converts wavelength (A) to theta angle (deg) """
  return 180./np.pi*np.arcsin(wavelength/2/gSi111dspacing )

def alioToWavelength(alio):
  """ Function that converts alio position (mm) to wavelength (A) """
  theta = alioToTheta(alio)
  return thetaToWavelength(theta)

def alioToE(alio):
  """ Function that converts alio position (mm) to photon energy (keV) """
  l = alioToWavelength(alio)
  return wavelengthToE(l)

def wavelengthToE(wavelength):
  """ Fucntion that converts wavelength (A) to photon energy (keV) """
  return 12.39842/wavelength

def EToWavelength(E):
  """ Function that converts photon energy (keV) to wavelength (A) """
  return 12.39842/E

def EToAlio(E):
  """ Function that converts photon energy (keV) to alio position (mm) """
  l = EToWavelength(E)
  t = wavelengthToTheta(l)
  alio = thetaToAlio(t)
  return alio

