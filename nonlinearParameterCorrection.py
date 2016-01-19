""" This module is to correct CSPAD images as discussed
in J Synchrotron Radiat. 2015 May 1; 22(Pt 3): 584–591.
doi:  10.1107/S1600577515005536
Correction of complex nonlinear signal response from a pixel array detector"""
from ixppy import tools,wrapFunc,dataset,matchEvents
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
  

#################3
def getCorrectionFunc(dmat=None,i=None,ic=None,order=5,sc=None,search_dc_limits=None,corrtype='corrNonLin',wrapit=True):
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
  if ic is None:
    ic = np.mean(i)
  if search_dc_limits is not None:
    search_dc_limits = iterfy(search_dc_limits)
    if len(search_dc_limits)==1:
      msk = (i>i-np.abs(search_dc_limits)) & (i<i+np.abs(search_dc_limits))
    elif len(search_dc_limits)==2:
      msk = (i>i-np.min(search_dc_limits)) & (i<i+np.max(search_dc_limits))
    p0 = tools.polyFit(i[msk],dmat[msk,...],2)

    dc = tools.polyVal(p0,i0_wp)
    pc = tools.polyFit(i-ic,Imat-dc,order,removeOrders=[0])
    if corrtype is 'corrNonLin':
      pcprime = tools.polyDer(pc)
    #c = lambda(i): polyVal(pc,i-ic) + dc
    def c(i):
      #print np.shape(pc)
      return tools.polyVal(pc,i-ic) + dc
  else:
    pc = tools.polyFit(i-ic,dmat,order,removeOrders=[])
    #c = lambda(i): tools.polyVal(pc,i-ic)
    def c(i):
      #print np.shape(pc)
      return tools.polyVal(pc,i-ic)
    dc = c(ic)
    if corrtype is 'corrNonLin':
      pcprime = tools.polyDer(pc)
  if corrtype is 'corrNonLin':
    pcprime = tools.polyDer(pc)
    c_prime = lambda(i): tools.polyVal(pcprime,i-ic)
    cprimeic = c_prime(ic)
    if sc is None:
      sc = c(ic)
    def corrFunc(Dm,im):
      im = np.asarray(im).ravel()
      return (sc.T/ic*im).T + cprimeic/c_prime(im)* sc/dc * (Dm-c(im)) #!!!
  elif corrtype is 'removeDep':
    def corrFunc(Dm,im):
      im = np.asarray(im).ravel()
      return Dm-c(im)+dc


  #return corrFunc
  if wrapit:
    def corrFuncTransposed(Dm,im=None,normalize=False,fillValue=np.nan):
      if im is None:
	im = np.apply_over_axes(np.nansum,Dm,range(1,np.ndim(Dm))).ravel()
      cr = corrFunc(Dm,im)
      if normalize:
	im = im.ravel()
	for dimno in range(cr.ndim-1):
	  im = np.expand_dims(im,-1)
	cr/=im
      cr[~np.logical_and(im>np.min(i),im<np.max(i)),...] *= fillValue
      return cr
    #else: 
	#return corrFunc(D.swapaxes(0,-1),i).swapaxes(0,-1)

    corrFuncWrapped = wrapFunc(corrFuncTransposed,transposeStack=False)
    def corrFuncWrap(Dm,im=None,normalize=False,fillValue=np.nan):
      if im is not None:
	Df = Dm*im.filter([np.min(i),np.max(i)]).ones()
      else:
	Df = Dm
      return corrFuncWrapped(Df,im=im,normalize=normalize,fillValue=fillValue)
    return corrFuncWrap
  else:
    return corrFunc


################3333
def getCorrectionFunc_old(order=5,Imat=None,i0=None,i0_wp=None,fraclims_dc=[.9,1.1],wrapit=True):
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
  if i0_wp is None:
    i0_wp = np.mean(i0)
  msk = tools.filtvec(i0,i0_wp*np.asarray(fraclims_dc))
  p0 = tools.polyFit(i0[msk],Imat[msk,...],2)
  dc = tools.polyVal(p0,i0_wp)
  comps = tools.polyFit(i0-i0_wp,Imat-dc,order,removeOrders=[0])
  compsder = tools.polyDer(comps)
  c = lambda(i): tools.polyVal(comps,i-np.asarray(tools.iterfy(i0_wp)))+dc
  c_prime = lambda(i): tools.polyVal(compsder,i-np.asarray(tools.iterfy(i0_wp)))
  t = lambda(i): (c_prime(i0_wp).T * (i-i0_wp)).T + dc
  
  cprimeic = c_prime(i0_wp)
  dcorr_const = -cprimeic*i0_wp + c(i0_wp) - t(0) 
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




class CorrNonLin(object):
  def __init__(self,data=None,I0=None,Imat=None,Iref=None,fina=None):
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
    if Iref is not None:
      self.Iref,self.data = matchEvents(Iref,data)
    else:
      self.data = data

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
    fina = 'tmp_getRefIntensity_' \
	+ datetime.datetime.now().isoformat() + '.ixp.h5'
    print fina
    self.Iref.setFile(fina)
    self.Iref.evaluate()
    self.Iref = self.Iref.get_memdata()[0]
    os.remove(fina)

  def getI0Imat(self,bins=None,evaluate=False):
    digi = (self.refDataFilter*self.Iref).digitize(bins=bins)
    self.Iref_good = self.Iref * digi.ones()
    self.I0 = digi.scan.bincenters
    self.Imat = 1/digi*self.data*self.I0
    if evaluate:
      fina = 'tmp_getImat_' \
	  + datetime.datetime.now().isoformat() + '.ixp.h5'
      print fina
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
  def getCorrFunc(self,order=5,ic=None,search_dc_limits=None, wrapit=True):
    if not ((self.Imat is None) and (self.I0 is None)):
      self.correct = getCorrectionFunc(order=order,dmat=self.Imat,i=self.I0,ic=ic,search_dc_limits=search_dc_limits, wrapit=wrapit, corrtype='corrNonLin')
      return self.correct

  def testCorrfunc(self,order=5,ic=None):
    fig = tools.nfigure('test_correction_func_order_%d'%order)
    plt.clf()
    #fig,ax = plt.subplots(1,2,num=fig.number)
    fig,ax = plt.subplots(1,2,num=fig.number,sharex=True,sharey=True)
    plt.axes(ax[0])
    it = (self.Imat.T/self.I0).T
    tools.imagesc(np.arange(np.prod(np.shape(self.Imat)[1:])),self.I0,((it/np.mean(it,0))-1).reshape(it.shape[0],-1))
    tools.clim_std(2)
    plt.colorbar()
    plt.draw()
    cf = self.getCorrFunc(order=order,ic=ic,wrapit=False)
    Icorr = cf(self.Imat,self.I0)
    plt.axes(ax[1])
    it = (Icorr.T/self.I0).T
    tools.imagesc(np.arange(np.prod(np.shape(self.Imat)[1:])),self.I0,((it/np.mean(it,0))-1).reshape(it.shape[0],-1))
    tools.clim_std(2)
    plt.colorbar()



class CorrNonLinDep(object):
  def __init__(self,data=None,Iref=None,I0=None,Imat=None,fina=None):
    self.data = data
    self.refDataFilter = 1
    self.Iref = Iref
    if fina is not None:
      assert fina[-7:]=='.ixp.h5', "File name has to be of extension ... .ixp.h5"
      self.dataset = dataset(fina)
      if 'corrNonLinRem_I0' in self.dataset.__dict__:
	self.I0 = self.dataset['corrNonLinRem_I0']
      else:
	self.I0 = None
      if 'corrNonLinRem_Imat' in self.dataset.__dict__:
	self.Imat = self.dataset['corrNonLinRem_Imat']
      else:
	self.Imat = None
    else:
      self.dataset = None


  def getRefdataMask(self,*args):
    flt = 1
    if 'step' in args:
      flt *= (self.data.ones()*self.data.scan[0]).filter().ones()
    fig,ax = plt.subplots(1,2,num=fig.number,sharex=True,sharey=True)
    tools.imagesc(((it/np.mean(it,0))-1).reshape(it.shape[0],-1))
    self.refDataFilter = flt

  def _getRefdata(self):
    return self.refDataFilter*self.data
  refData = property(_getRefdata)

  def getRefIntensity(self,imagemask=None):
    self.Iref = self.refDataFilter * nansum(self.data)
    fina = 'tmp_getRefIntensity_' \
	+ datetime.datetime.now().isoformat() + '.ixp.h5'
    print fina
    self.Iref.setFile(fina)
    self.Iref.evaluate()
    self.Iref = self.Iref.get_memdata()[0]
    os.remove(fina)

  def getI0Imat(self,bins=None,evaluate=False):
    digi = (self.refDataFilter*self.Iref*self.data.ones()).digitize(bins=bins)
    self.I0 = digi.scan.bincenters
    self.Imat = digi.ones()*self.data
    if evaluate:
      fina = 'tmp_getImat_' \
	  + datetime.datetime.now().isoformat() + '.ixp.h5'
      print fina
      self.Imat.setFile(fina)
      self.Imat.evaluate()
      self.Imat = np.asarray(self.Imat.median())
      os.remove(fina)
    else:
      self.Imat = np.asarray(self.Imat.median())
      if self.dataset is not None:
	self.dataset['corrNonLinRem_Imat'] = self.Imat
	self.dataset['corrNonLinRem_I0'] = self.I0
	self.dataset.save()

  def getCorrFunc(self,order=5,ic=None,search_dc_limits=None, wrapit=True):
    if not ((self.Imat is None) and (self.I0 is None)):
      self.correct = getCorrectionFunc(order=order,dmat=self.Imat,i=self.I0,ic=ic,search_dc_limits=search_dc_limits, wrapit=wrapit, corrtype='removeDep' )
      return self.correct

  def testCorrfunc(self,order=5,ic=None):
    fig = tools.nfigure('test_correction_func_order_%d'%order)
    plt.clf()
    fig,ax = plt.subplots(1,2,num=fig.number,sharex=True,sharey=True)
    plt.axes(ax[0])
    it = (self.Imat.T/self.I0).T
    tools.imagesc(np.arange(np.prod(np.shape(self.Imat)[1:])),self.I0,((it/np.mean(it,0))-1).reshape(it.shape[0],-1))
    #tools.imagesc(np.arange(np.shape(self.Imat)[1]),self.I0,(it/np.mean(it,0))-1)
    tools.clim_std(2)
    plt.colorbar()
    plt.draw()
    cf = self.getCorrFunc(order=order,ic=ic,wrapit=False)
    Icorr = cf(self.Imat,self.I0)
    plt.axes(ax[1])
    it = Icorr
    #tools.imagesc((it/np.mean(it,0))-1)
    tools.imagesc(np.arange(np.prod(np.shape(self.Imat)[1:])),self.I0,((it/np.mean(it,0))-1).reshape(it.shape[0],-1))
    tools.clim_std(2)
    plt.colorbar()

