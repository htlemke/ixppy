import numpy as np
import scipy as sp
import numpy.linalg as linalg
from scipy.special import erf,gamma
from scipy.stats import skew
from scipy.interpolate import interp1d 
from scipy import optimize
import matplotlib
import pylab as pl
import time
import types
import numpy.ma as ma
import toolsDistrAndHist
from toolsVarious import iterfy
try:
  import minuit
except:
  print "Minuit not found"



def gaussConvExp(x,x0,sig,tau,amp,y0):
  fval = amp*.5*(np.exp(((np.sqrt(2)*sig)**2-4*(x-x0)*tau)/(4*tau**2)) *  erf((np.sqrt(2)*sig)/(2*tau)-(x-x0)/(np.sqrt(2)*sig)) + erf((x-x0)/(np.sqrt(2)*   sig)) - np.exp(((np.sqrt(2)*sig)**2-4*(x-x0)*tau)/(4*tau**2))+1)+y0
  return fval

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



def gauss_par(par=dict(A=[],pos=[],sig=[]),dat=[]):
	res = par['A']*np.exp(-(dat-par['pos'])**2/par['sig']**2/2)
	return res

def gauss(dat,A,pos,sig):
	res = A*np.exp(-(dat-pos)**2/sig**2/2)
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
  g = toolsDistrAndHist.gauss_norm
  o = np.empty_like(y)
  from scipy.integrate import simps
  for i in range(len(x)):
    w = g(par={"sig": sig, "pos": x[i], "area": 1},xdat=x)
    o[i] = simps(w*y,x)/simps(w,x)
  return o

def pearsonVIIsplit(xdata,Amp,xc,H,A,mL,mH):
  #% Usage: [ycalc] = pearsonVIIsplit(X, xdata)
  #% Split pearson VII peak profile as described in 
  #% Toraya, H. (1990). "Array-Type Universal Profile 
  #% Function for Powder Pattern Fitting." 
  #% Journal of Applied Crystallography 23: 485-491. 
  #%
  #% I would be interested in another version where the profile is
  #% ANALYTICALLY convoluted with a rectangular function of certain width.
  #% I tried one day, it should work, but I gave up.
  #% 
  #% Henrik T. Lemke, March 2009
  # Amp=height
  # xc = center
  # H = FWHM
  # A = assymetry, symmetric is 1
  # mL,mH = wings >=1

  xdata = np.asarray(iterfy(xdata))

  ycalc = np.zeros_like(xdata)
  ycalc[xdata<=xc] = Amp*2.*(1+A)/(np.sqrt(np.pi)*H)\
      *((A*gamma(mL-.5))/(np.sqrt(2.**(1./mL)-1)*gamma(mL)) \
	  + gamma(mH-.5)/(np.sqrt(2.**(1./mH)-1)*gamma(mH)))**(-1)\
	  *(1 + (2.**(1./mL)-1)*((1+A)/A)**2 * ((xdata[xdata<=xc]-xc)/H)**2) **(-mL)
  A = 1/A
  ycalc[xdata>xc] = Amp*2.*(1+A)/(np.sqrt(np.pi)*H)\
      *((A*gamma(mH-.5))/(np.sqrt(2.**(1./mH)-1)*gamma(mH)) \
	  + gamma(mL-.5)/(np.sqrt(2.**(1./mL)-1)*gamma(mL)))**(-1)\
	  *(1 + (2.**(1./mH)-1)*((1+A)/A)**2 * ((xdata[xdata>xc]-xc)/H)**2) **(-mH)
  return ycalc

def peakAna(x,y,nb=3,plotpoints=False):
	""" nb = number of point (on each side) to use as background"""
	## get background
	xb = np.hstack((x[0:nb],x[-(nb):]))
	yb = np.hstack((y[0:nb],y[-(nb):]))
	a = np.polyfit(xb,yb,1)
	b = np.polyval(a,x)
	yf = y-b
	yd = np.diff(yf)

	## determine whether peak or step
	ispeak = np.abs(skew(yf))>np.abs(skew(yd))
	if ispeak:
		yw = yf
		xw = x
	else:
		yw = yd
		xw = (x[1:]+x[0:-1])/2
		## get background
		xwb = np.hstack((xw[0:nb],xw[-(nb):]))
		ywb = np.hstack((yw[0:nb],yw[-(nb):]))
		aw = np.polyfit(xwb,ywb,1)
		bw = np.polyval(aw,xw)
		yw = yw-bw
	

	Iw = (xw[1:]-xw[0:-1])*(yw[1:]+yw[0:-1])/2
	if sum(Iw)<0:
		yw = -yw

	## get parameters	
	mm = yw.argmax(0)
	PEAK = xw[mm]
	ywmax = yw[mm]
	gg = (yw[:mm][::-1]<(ywmax/2)).argmax()
	ip = interp1d(yw.take([mm-gg-1,mm-gg]),xw.take([mm-gg-1,mm-gg]),kind='linear')
	xhm1 = ip(ywmax/2)
	gg = (yw[mm:]<(ywmax/2)).argmax()
	ip = interp1d(yw.take([mm+gg,mm+gg-1]),xw.take([mm+gg,mm+gg-1]),kind='linear')
	xhm2 = ip(ywmax/2)

	FWHM = np.abs(xhm2-xhm1)
	CEN = (xhm2+xhm1)/2
	if plotpoints and ispeak is True:
		# plot the found points for center and FWHM edges
		ion()
		pl.hold(True)
		pl.plot(x,b,'g--')
		pl.plot(x,b+ywmax,'g--')

		pl.plot([xhm1,xhm1],polyval(a,xhm1)+[0,ywmax],'g--')
		pl.plot([xhm2,xhm2],polyval(a,xhm2)+[0,ywmax],'g--')
		pl.plot([CEN,CEN],polyval(a,CEN)+[0,ywmax],'g--')
		pl.plot([xhm1,xhm2],[polyval(a,xhm1),polyval(a,xhm2)]+ywmax/2,'gx')


		pl.draw()

	
	if not ispeak:
	  try:
	        # findings start of step coming from left.

		std0 = sp.std(y[0:nb])
		nt = nb
	
		while (sp.std(y[0:nt])<(2*std0)) and (nt<len(y)):
			nt = nt+1
		lev0 = sp.mean(y[0:nt])
	
		# findings start of step coming from right.
		std0 = sp.std(y[-nb:])
		nt = nb

		while (sp.std(y[-nt:])<(2*std0)) and (nt<len(y)):
			nt = nt+1
		lev1 = sp.mean(y[-nt:])
                

		gg = np.abs(y-((lev0+lev1)/2)).argmin()     

		ftx = y[gg-2:gg+2]
		fty  = x[gg-2:gg+2]
		if ftx[-1]<ftx[0]:
			ftx = ftx[::-1]
			fty = fty[::-1]
	
		ip = interp1d(ftx,fty,kind='linear')
		CEN = ip((lev0+lev1)/2)
		
		gg = np.abs(y-(lev1+(lev0-lev1)*0.1195)).argmin()     

		ftx = y[gg-2:gg+2]
		fty  = x[gg-2:gg+2]
		if ftx[-1]<ftx[0]:
			ftx = ftx[::-1]
			fty = fty[::-1]
		#print " %f %f %f %f %f" % (ftx[0],ftx[1],fty[0],fty[1],lev1+(lev0-lev1)*0.1195)
		ip = interp1d(ftx,fty,kind='linear')
		H1 = ip((lev1+(lev0-lev1)*0.1195))
		#print "H1=%f" % H1

		gg = np.abs(y-(lev0+(lev1-lev0)*0.1195)).argmin()     

		ftx = y[gg-2:gg+2]
		fty  = x[gg-2:gg+2]
		
		if ftx[-1]<ftx[0]:
			ftx = ftx[::-1]
			fty = fty[::-1]
#		print " %f %f %f %f %f" % (ftx[0],ftx[1],fty[0],fty[1],lev0+(lev1-lev0)*0.1195)
		ip = interp1d(ftx,fty,kind='linear')
		H2 = ip((lev0+(lev1-lev0)*0.1195))
		#print "H2=%f" % abs(H2-H1)
		FWHM = abs(H2-H1)
		if plotpoints is True:
			# plot the found points for center and FWHM edges
			pl.ion()
			pl.hold(True)
			pl.plot([x.min(),x.max()],[lev0,lev0],'g--')
			pl.plot([x.min(),x.max()],[lev1,lev1],'g--')

			pl.plot([H2,H2],[lev0,lev1],'g--')
			pl.plot([H1,H1],[lev0,lev1],'g--')
			pl.plot([CEN,CEN],[lev0,lev1],'g--')
			pl.plot([H2,CEN,H1],[lev0+(lev1-lev0)*0.1195,(lev1+lev0)/2,lev1+(lev0-lev1)*0.1195],'gx')


			pl.draw()
	  except:
	        CEN = np.nan
	        FWHM = np.nan
	        PEAK = np.nan

	return (CEN,FWHM,PEAK)

def fitCircle(x,y,w=1.):
  def calc_R(x,y, xc, yc):
      """ calculate the distance of each 2D points from the center (xc, yc) """
      return np.sqrt((x-xc)**2 + (y-yc)**2)
   
  def f(c, x, y, w):
      """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
      Ri = calc_R(x, y, *c)
      return w*(Ri - Ri.mean())
   
  x_m = np.mean(x)
  y_m = np.mean(y)
  center_estimate = x_m, y_m
  center, ier = optimize.leastsq(f, center_estimate, args=(x,y,w))
  xc, yc = center
  Ri       = calc_R(x, y, *center)
  R        = Ri.mean()
  residu   = np.sum((Ri - R)**2)
  return xc, yc, R, residu
