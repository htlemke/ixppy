import numpy as np
import time
import types
import numpy.ma as ma

def eV2reccm(eV):
	reccm =  eV* 8065.54445
	return reccm

def reccm2eV(reccm):
	eV = reccm / 8065.54445
	return eV

def eV2nm(eVvec):
	nmvec = 1e9*h_planck()* c_light() / eV2J(eVvec)
	return nmvec

def nm2eV(nmvec):
	eVvec = J2eV(1e9*h_planck()*c_light()/nmvec)
	return eVvec

def eV2J(eV):
	J = 1.60217646e-19 * eV
	return J

def J2eV(J):
	eV = J/1.60217646e-19
	return eV

def c_light():
	c = 299792458 # m/s
	return c

def h_planck():
	h = 6.626068e-34 # m2 kg / s
	return h

def E2lam(E):
	lam = 12.39842 /E
	return lam

def lam2E(lam):
	E = 12.39842 / lam;  #/keV
	return E
