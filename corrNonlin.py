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

def corrNonlin(expar,Iuncorr,comps,exparWP,Iwp,):
  pol = np.vander(expar-exparWP,order+1)
  
  if Iwp is not None:
    np.delete(pol,-1,1)
  


def removeNonlin(components,i0s):
  pol = np.vander(i0,np.shape(components)[0])
  pol = np.hstack([pol[:,:-2],pol[:,-1]])
  return np.array(np.matrix(pol)*np.matrix(components))

def getCorr(order=5,i0=None,Imat=None,i0_wp=1e6,fraclims_dc=[.9,1.1]):
  """ Getting nonlinear correction factors form a calibration dataset consiting of:
  i0     	array of intensities the calibration has been made for
  Imat   	2D array of the corresponding reference patterns, in each row 
  		there is one ravelled array of each intensity bin in i0.
  i0_wp		a working point around which a correction polynomial will be
  		developed for each pixel.
  order		the polynomial order up to which will be deveoped.
  fraclims_dc	relative factor for the i0,Imat data limits which are used to
  		determine the working point location.
  
  Returns

  """

  #i0,Imat = getData()
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
  def dcorr(i,D):
    return (i*cprimeic.T + dcorr_const.T + ((D-c(i))*cprimeic/c_prime(i)).T).T
    #return (i*cprimeic.T + dcorr_const.T ).T
  return dcorr,comps,t





  tools.nfigure('testplot')
  plt.clf()
  plt.subplot(1,2,1)
  Imean = (Imat.T/i0).T
  tools.imagesc(np.asarray([ti / np.mean(Imean[-10:,:],0) for ti in Imean]))
  tools.clim_std(6)
  cl = plt.gci().get_clim()
  plt.colorbar()
  plt.set_cmap(plt.cm.RdBu_r)
  
  
  plt.subplot(1,2,2)

  cmps = copy.copy(comps)
  cmps[-2,:] = 0
  cc = lambda(i): tools.polyVal(cmps,i-np.asarray(tools.iterfy(i0_wp)))
  Ir = Imat-c(i0)+t(i0)-t(0)
  Ir = dcorr(i0,Imat)
  #Ir = ((Imat-cc(i0)).T/i0).T
  #tools.imagesc(Ir) 
  Ir = (Ir.T/i0).T
  tools.imagesc(np.asarray([ti / np.mean(Ir[-10:,:],0) for ti in Ir]))
  plt.clim(cl)
  plt.colorbar()
  plt.set_cmap(plt.cm.RdBu_r)
  plt.draw()

  tools.nfigure('testplot_components')
  plt.clf()
  ah = None
  for n,comp in enumerate(comps):
    if ah is None:
      ah = plt.subplot(len(comps),1,n+1)
    else:
      plt.subplot(len(comps),1,n+1,sharex=ah)
    plt.plot(comp)
    lims = np.percentile(comp,[1,99])

    plt.ylim(lims)


  return c,c_prime 
  #dci = 

def rearrangeData(data,mask):
  re = np.ones((len(data),len(mask)))*np.nan
  msk = mask.nonzero()
  for n,dat in enumerate(data):
    re[n].put(msk,dat)
  return re.reshape([len(data),32,185,388])

def makeplots(i0,Imat,mask,ind=[27, 29, 37, 67],order=6):
  Imatref = Imat[np.ix_([27, 29, 37, 67])]
  i0ref = i0[np.ix_([27, 29, 37, 67])]
  dcorr,comps,t = getCorr(order,i0,Imat,i0ref[0])
  patt = cspad.CspadPattern(Nx=300,Ny=300)
  
  if 1:
    Dcorr = dcorr(i0ref,Imatref)
    DcorrNorm = (Dcorr.T/i0ref).T
    DcorrNorm = rearrangeData(DcorrNorm,mask)
    ImatrefNorm = (Imatref.T/i0ref).T
    ImatrefNorm = rearrangeData(ImatrefNorm,mask)
    fig = tools.nfigure('figure2')
    #fig.set_size_inches(3.5,3)
    plt.clf()
    ah = None
    for n in range(3):
      if ah is None:
        ah = plt.subplot(2,3,n+1)
	tah = ah
      else:
        tah = plt.subplot(2,3,n+1,sharex=ah,sharey=ah)
      patt.imageShow((ImatrefNorm[1+n]/ImatrefNorm[0]-1)*100)
      plt.set_cmap(plt.cm.RdBu_r)
      plt.clim([-5,5])
      plt.axis('equal')
      plt.axis([60000,13e4,6e4,13e4])
      plt.setp(tah.get_xticklabels(), visible=False) 
      plt.setp(tah.get_yticklabels(), visible=False) 
      plt.text(6.4e4,11.8e4,['(a)','(b)','(c)'][n],fontsize=20,bbox={'facecolor':'white', 'alpha':1, 'pad':10})


      tah = plt.subplot(2,3,n+4,sharex=ah,sharey=ah)
      pp = patt.bin((DcorrNorm[1+n]/DcorrNorm[0]-1)*100)
      im = tools.imagesc(patt.xVec,patt.yVec,pp)
      plt.set_cmap(plt.cm.RdBu_r)
      plt.clim([-5,5])
      plt.axis('equal')
      plt.axis([60000,13e4,6e4,13e4])
      plt.setp(tah.get_xticklabels(), visible=False) 
      plt.setp(tah.get_yticklabels(), visible=False) 
      plt.text(6.4e4,11.8e4,['(d)','(e)','(f)'][n],fontsize=20,bbox={'facecolor':'white', 'alpha':1, 'pad':10})

    fig.subplots_adjust(left=.05,top=.95,bottom=.05,right=0.85,hspace=.05,wspace=.05)
    cbar_ax = fig.add_axes([0.88, 0.05, 0.03, 0.9])
    fig.colorbar(im, cax=cbar_ax,label='Percent')

  # FIGURE 4

  fig = tools.nfigure('figure4')
  #fig.set_size_inches(3.5,1.5)
  plt.clf()
  ah = None
  count=1
  toplot = [t(0),comps[-3],comps[-4]]
  for n,comp in enumerate(toplot):
    if ah is None:
      ah = plt.subplot(1,3,count)
      tah = ah
    else:
      tah = plt.subplot(1,3,count,sharex=ah,sharey=ah)
    timg = rearrangeData([comp],mask)
    patt.imageShow(timg)
    plt.set_cmap(plt.cm.RdBu_r)
    #plt.clim([-5,5])
    plt.axis('equal')
    plt.axis([60000,13e4,6e4,13e4])
    plt.setp(tah.get_xticklabels(), visible=False) 
    plt.setp(tah.get_yticklabels(), visible=False) 
    plt.text(6.4e4,11.8e4,['$t(0)$','$g=2$','$g=3$','d','e','f','g','h','i','j','k'][n],fontsize=20,bbox={'facecolor':'white', 'alpha':1, 'pad':10})
    count += 1

  fig.subplots_adjust(left=.05,top=.95,bottom=.05,right=0.95,hspace=.05,wspace=.05)


  #fig.subplots_adjust(right=0.85,hspace=.05,wspace=.05)
  #cbar_ax = fig.add_axes([0.9, 0.1, 0.02, 0.8])
  #fig.colorbar(im, cax=cbar_ax,ylabel=')

