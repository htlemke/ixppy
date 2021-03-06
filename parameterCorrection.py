from toolsLog import logbook
from ixppy import tools,wrapFunc,dataset
import pylab as plt
import numpy as np
from scipy import linalg,io
import copy
from toolsExternalWrapped import nansum
import datetime,os



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
  

def getCorrectionFunc(dmat=None,i=None,ic=None,order=5,sc=None,search_dc_limits=None, removeDependence=True):
  """ 
  Create nonlinear correction function from a calibration dataset consiting of:
    i     	array of intensity values (floats) of the calibration
    dmat   	ND array of the reference patterns corresponding to values of i,
                The first dimension corresponds to the calibration intensity 
		values and has the same length as i.
    ic		A working point around which a polynomial correction will be
  		developed for each pixel.
    order	the polynomial order up to which the correction will be 
                deveoped.
    sc		optional: calibration image at ic. default is the image at ic
    search_dc_limits	absolute limits around ic which are used to determine the 
    		calibration value of ic as linear approximation of a short interval. 
		optional, can sometimes help to avoid strong deviations of the 
		polynomial approximatiuon from the real measured points.
  
  Returns corrFunc(D,i), a function that takes an ND array input for correction
                (1st dimension corresponds to the different intensity values) 
		as well as the intensity array i.
  """
  if search_dc_limits is not None:
    search_dc_limits = iterfy(search_dc_limits)
    if len(search_dc_limits)==1:
      msk = (i>i-np.abs(search_dc_limits)) & (i<i+np.abs(search_dc_limits))
    elif len(search_dc_limits)==2:
      msk = (i>i-np.min(search_dc_limits)) & (i<i+np.max(search_dc_limits))
    p0 = tools.polyFit(i[msk],dmat[msk,...],2)
    dc = tools.polyVal(p0,i0_wp)
    pc = tools.polyFit(i-ic,Imat-dc,order,removeOrders=[0])
    pcprime = tools.polyDer(pc)
    c = lambda(i): polyVal(pc,i-ic) + dc
  else:
    pc = polyFit(i-ic,dmat,order,removeOrders=[])
    pcprime = polyDer(pc)
    c = lambda(i): polyVal(pc,i-ic)
    dc = c(ic)
  c_prime = lambda(i): polyVal(pcprime,i-ic)
  cprimeic = c_prime(ic)
  if sc is None:
    sc = c(ic)
  def corrFunc(D,i):
    i = np.asarray(i).ravel()
    return (sc.swapaxes(0,-1)/ic*i).swapaxes(0,-1) + cprimeic/c_prime(i)* sc/dc * (D-c(i))
  def remDepFunc(D,i):
    i = np.asarray(i).ravel()
    return D-c(i)+dc
  if removeDependence:
    corrFunc = remDepFunc

  if wrapit:
    def corrFuncTransposed(D,i,normalize=False,fillValue=np.nan):
      cr = corrFunc(D.swapaxes(0,-1),i).swapaxes(0,-1)
      if normalize:
	cr/=i
      cr[:,~np.logical_and(i>np.min(i0),i<np.max(i0))] *= fillValue
      return cr

    corrFuncWrapped = wrapFunc(corrFuncTransposed,transposeStack=True)
    def corrFuncWrap(D,i,normalize=False,fillValue=np.nan):
      Df = D*i.filter([np.min(i0),np.max(i0)]).ones()
      return corrFuncWrapped(Df,i=i,normalize=normalize,fillValue=fillValue)
    return corrFuncWrap
  else:
    return corrFunc

  return corrFunc

def getCorrectionFunc(order=5,Imat=None,p=None,pc=None,fraclims_dc=[.9,1.1],wrapit=True):
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
  if pc is None:
    pc = np.mean(p)
  msk = tools.filtvec(p,pc*np.asarray(fraclims_dc))
  p0 = tools.polyFit(p[msk],Imat[msk,...],2)
  dc = tools.polyVal(p0,pc)
  comps = tools.polyFit(p-pc,Imat-dc,order,removeOrders=[0])
  compsder = tools.polyDer(comps)
  c = lambda(i): tools.polyVal(comps,i-np.asarray(tools.iterfy(pc)))+dc
  c_prime = lambda(i): tools.polyVal(compsder,i-np.asarray(tools.iterfy(pc)))
  t = lambda(i): (c_prime(pc).T * (i-pc)).T + dc
  
  cprimeic = c_prime(pc)
  dcorr_const = -cprimeic*pc + c(pc) - t(0) 
  def corrFunc(D,i):
    i = i.ravel()
    return cprimeic * ( i + ((D-c(i))/c_prime(i)).swapaxes(0,-1) ).swapaxes(0,-1)

  if wrapit:
    def corrFuncTransposed(D,i=None,normalize=False,fillValue=np.nan):
      if i is None:
	i = np.apply_over_axes(np.nansum,D,range(np.ndim(D)-1)).ravel()
      cr = corrFunc(D.swapaxes(0,-1),i).swapaxes(0,-1)
      if normalize:
	cr/=i
      cr[:,~np.logical_and(i>np.min(i0),i<np.max(i0))] *= fillValue
      return cr
    #else: 
	#return corrFunc(D.swapaxes(0,-1),i).swapaxes(0,-1)

    corrFuncWrapped = wrapFunc(corrFuncTransposed,transposeStack=True)
    def corrFuncWrap(D,i=None,normalize=False,fillValue=np.nan):
      if i is not None:
	Df = D*i.filter([np.min(i0),np.max(i0)]).ones()
      else:
	Df = D
      return corrFuncWrapped(Df,i=i,normalize=normalize,fillValue=fillValue)
    return corrFuncWrap
  else:
    return corrFunc
  #return corrFunc




class CorrPar(object):
  def __init__(self,data=None,Iref=None,Imat=None,fina=None):
    self.data = data
    self.Iref = Iref
    self.refDataFilter = 1
    if fina is not None:
      assert fina[-7:]=='.ixp.h5', "File name has to be of extension ... .ixp.h5"
      self.dataset = dataset(fina)
      if 'corrNonLin_I0' in self.dataset.__dict__:
	self.I0 = self.dataset['corrNonLin_I0']
      else:
	self.I0 = None
      if 'corrNonLin_Imat' in self.dataset.__dict__:
	self.Imat = self.dataset['corrNonLin_Imat']
      else:
	self.Imat = None
    else:
      self.dataset = None

  def getRefdataMask(self,*args):
    flt = 1
    if 'step' in args:
      flt *= (self.data.ones()*self.data.scan[0]).filter().ones()
    self.refDataFilter = flt

  def _getRefdata(self):
    return self.refDataFilter*self.data
  refData = property(_getRefdata)

  def getRefIntensity(self,imagemask=None):
    self.Iref = self.refDataFilter * nansum(self.data)
    fina = 'tmp_getRefIntensity_' + \
	  datetime.datetime.now().isoformat() + '.ixp.h5'
    logbook(fina)
    self.Iref.setFile(fina)
    self.Iref.evaluate()
    self.Iref = self.Iref.get_memdata()[0]
    os.remove(fina)

  def getI0Imat(self,bins=None,evaluate=False):
    digi = (self.refDataFilter*self.Iref).digitize(bins=bins)
    self.I0 = digi.scan.bincenters
    self.Imat = digi.ones()*self.data
    if evaluate:
      fina = 'tmp_getImat_' \
	  + datetime.datetime.now().isoformat() + '.ixp.h5'
      logbook(fina)
      self.Imat.setFile(fina)
      self.Imat.evaluate()
      self.Imat = np.asarray(self.Imat.mean())
      os.remove(fina)
    else:
      self.Imat = np.asarray(self.Imat.mean())
      if self.dataset is not None:
	self.dataset['corrNonLin_Imat'] = self.Imat
	self.dataset['corrNonLin_I0'] = self.I0
	self.dataset.save()
  def getCorrFunc(self,order=5,i0_wp=None,fraclims_dc=[.9,1.1], wrapit=True):
    if not ((self.Imat is None) and (self.I0 is None)):
      self.correct = getCorrectionFunc(order=order,Imat=self.Imat,i0=self.I0,i0_wp=i0_wp,fraclims_dc=fraclims_dc, wrapit=wrapit)
      return self.correct

  def testCorrfunc(self,order=5,ic=None):
    fig = tools.nfigure('test_correction_func_order_%d'%order)
    plt.clf()
    fig,ax = plt.subplots(1,2,num=fig.number)
    plt.axes(ax[0])
    it = (self.Imat.T/self.I0).T
    tools.imagesc(np.arange(np.shape(self.Imat)[1]),self.I0,(it/np.mean(it,0))-1)
    tools.clim_std(2)
    plt.colorbar()
    plt.draw()
    cf = self.getCorrFunc(order=order,i0_wp=ic,wrapit=False)
    Icorr = cf(self.Imat,self.I0)
    plt.axes(ax[1])
    it = (Icorr.T/self.I0).T
    tools.imagesc((it/np.mean(it,0))-1)
    tools.clim_std(2)
    plt.colorbar()









