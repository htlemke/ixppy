import numpy as np
import scipy as sp
import numpy.linalg as linalg
from scipy.special import erf
import matplotlib
import pylab as pl
import time
import types
import numpy.ma as ma
import tools
try:
  import minuit
except:
  print "Minuit not found"


def polyfitZero(x,y,degree=1):
  x = np.array(x).ravel()
  y = np.array(y).ravel()
  #thanks to Mark Mikofski
  z = np.zeros([len(x),degree])

  for n in range(degree):
    z[:,n] = x**(degree-n)

  pZero,_,_,_ = np.linalg.lstsq(z,y.T)
  #pZero = [pZero;0]';

      #end
  pZero = list(pZero)
  pZero.append(0)
  return np.array(pZero)


def chisqwrap(X,fun,xdat,ydat,bg_order=[]):
  """
  Usage e.g. with scipy.optimize.fmin:
  fmin(chisqwrap,[1,1,1,0,0],args = (gauss_amp,x,y,1))
  """
  ycalc = fun(X[0:np.shape(xdat)[0]-(bg_order)-1],xdat)+np.polyval(X[np.shape(xdat)[0]-(bg_order)-1:],xdat)
  chisq = np.sum((ydat-ycalc)**2)
  return chisq


def fitsimple(function,startpars,xdat,ydat,fixedpar=None):
  pardef = function.func_defaults[0]
  #if startpars is not dict:
    #splist = startpars
    #startpars = dict()
    #for tpar in pardef.keys():
      #startpars[tpar] = splist[pardef.keys().index(tpar)]

  def minfunc(x0):
    splist = x0
    startpars = dict()
    for tpar in pardef.keys():
      startpars[tpar] = splist[pardef.keys().index(tpar)]
    chsq = np.sum((ydat-function(startpars,xdat))**2)
    return chsq
    

  minpar = s.optimize.fmin(minfunc,startpars)

  return minpar


def minuitfit(function,startpars,xdat,ydat,edat=None,fixedpar=None,stepszpar=None,limitspar=None):
	"""
	Wrapper function for minuit to do Chi squared based fits on 
	simple function of format:
	ydat_calculated = function(dict(par1=...,par2=...,...) , xdat)
	The functions should be defined something like that
		def function(par=dict(p1=..,p2=..,...),xdata=[]):
			return par[p1]*xdata+par[p2]
	Then you can use it this way:
	par = dict(a=1,t0=0.,sig=0.03,tau=0.2,c=0.)
	m=mf(mystep,par,x,y,edat=e)
	m.simplex()
	m.migrad()
	fit = mystep(m.values,x)
	Other minuit parameters like the initial stepsize, parameter 
	fixing and limits (not implemented yet) can be given as both
	dictionary and list.
	HT Lemke 2011
	"""
	if edat is None:
		edat = np.ones(np.shape(ydat))
	elif np.shape(edat)==(1,):
		edat = edat*np.ones(np.shape(ydat))

	# internal global variables
	g = globals()
	g['idat'] = xdat
	g['odat'] = ydat
	g['edat'] = edat
	g['function'] = function

	# make dict from parameters is necessary and get most 
	# important strings for the minuit variable workaround
	pardef = function.func_defaults[0]
	varstr = ''
	varstrstartval = ''
	varstrval = ''
	if not (isinstance(startpars,dict)):
		splist = startpars
		startpars = dict()
		for tpar in pardef.keys():
			startpars[tpar] = splist[pardef.keys().index(tpar)]

	for tpar in pardef.keys():
		varstrstartval += '%s=startpars[\'%s\'],'%(tpar,tpar)
		varstrval += '%s=%s,'%(tpar,tpar)
		varstr += '%s,'%(tpar)
	varstr = varstr[:-1]
	varstrval = varstrval[:-1]
	varstrstartval = varstrstartval[:-1]



	# make fixed string
	if fixedpar is not None:
		if type(fixedpar) is dict:
			fixlist = []
			for tpar in pardef.keys():
				if tpar in fixedpar.keys():
					fixlist.append(fixedpar[tpar])
				else:
					fixlist.append(False)
		elif type(fixedpar) is list:
			fixlist = fixedpar
		fixstr = ''
		for tpar in pardef.keys():
			fixstr += 'fix_%s=%s,'%(tpar,fixlist[pardef.keys().index(tpar)])

		fixstr = fixstr[:-1]
	# make stepsize string
	if stepszpar is not None:
		if type(stepszpar) is dict:
			stepszlist = []
			for tpar in pardef.keys():
				if tpar in stepszpar.keys():
					stepszlist.append(stepszpar[tpar])
				else:
					fixlist.append(None)
		elif type(stepszpar) is list:
			stepszlist = stepszpar
		stepszstr = ''
		for tpar in pardef.keys():
			if stepszlist[pardef.keys().index(tpar)] is not None:
				stepszstr += 'err_%s=%s,'%(tpar,stepszlist[pardef.keys().index(tpar)])

		stepszstr = stepszstr[:-1]
	# make chisq function
	csstr  = 'def chisq('+varstr+'): cs = sum(((odat-function(dict('+varstrval+'),idat))/edat)**2); return cs'
	exec(csstr)
	mstr = 'minuit.Minuit(chisq,'+varstrstartval
	if stepszpar is not None:
		mstr += ","+stepszstr
	if fixedpar is not None:
		mstr += ","+fixstr
	mstr += ')'
	m = eval(mstr)
	return m


def gauss(par=dict(A=[],pos=[],sig=[]),dat=[]):
	res = par['A']*np.exp(-(dat-par['pos'])**2/par['sig']**2/2)
	return res

def heaviside(dat):
  hs = np.zeros_like(dat)
  hs[dat<0] = 0
  hs[dat>=0] = 1
  return hs

def convGauss2Exp(par=dict(pos=0,sig=.2,amp1=-.3,tau1=2,amp2=1,tau2=.3,y0=0),dat=np.linspace(-5,10,300)):
  x = dat
  res = par['amp1']*.5\
          *(np.exp(((np.sqrt(2)*par['sig'])**2-4*(x-par['pos'])*par['tau1'])/(4*par['tau1']**2))
            *erf((np.sqrt(2)*par['sig'])/(2*par['tau1'])-(x-par['pos'])/(np.sqrt(2)*par['sig'])) 
            + erf((x-par['pos'])/(np.sqrt(2)*par['sig'])) 
            - np.exp(((np.sqrt(2)*par['sig'])**2-4*(x-par['pos'])*par['tau1'])/(4*par['tau1']**2))+1)\
          +par['amp2']*.5\
          *(np.exp(((np.sqrt(2)*par['sig'])**2-4*(x-par['pos'])*par['tau2'])/(4*par['tau2']**2))
            *erf((np.sqrt(2)*par['sig'])/(2*par['tau2'])-(x-par['pos'])/(np.sqrt(2)*par['sig'])) 
            + erf((x-par['pos'])/(np.sqrt(2)*par['sig'])) 
            - np.exp(((np.sqrt(2)*par['sig'])**2-4*(x-par['pos'])*par['tau2'])/(4*par['tau2']**2))+1)\
          +par['y0']
  return res


def convolveGauss(x,y,sig):
  g = tools.gauss_norm
  o = np.empty_like(y)
  from scipy.integrate import simps
  for i in range(len(x)):
    w = g(par={"sig": sig, "pos": x[i], "area": 1},xdat=x)
    o[i] = simps(w*y,x)/simps(w,x)
  return o
