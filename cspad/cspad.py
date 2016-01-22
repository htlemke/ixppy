import pylab as pl
import numpy as np
import os,sys
from ixppy import tools,wrapFunc,Ixp
import copy
from ixppy.tools import nfigure,filtvec,poiss_prob,gauss_norm,polygonmask
#from functools import partial

parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
alignmentdir = parentdir+'/cspad/alignment'
#sys.path.insert(0,alignmentdir)
import alignment
from alignment import CalibPars          as calp
#import CalibParsEvaluated as cpe
#import CSPadConfigPars    as ccp
from alignment.CSPADPixCoords import CSPADPixCoords
import numpy as np
from matplotlib import pyplot as plt
import progressbar as pb
from toolsExternalWrapped import nansum
import datetime

g_maskUnbonded = dict()

def getCsPadPixCoordinates(path_calib=alignmentdir+'/calib-xpp-2013-01-29', 
                           rotation=0, 
                           mirror=False):
  if path_calib==None or path_calib=='newest':
    allfiles = os.listdir(alignmentdir)
    files = []
    for tf in allfiles:
      if 'calib-' in tf:
        files.append(tf)
    files.sort()

    if path_calib==None:
      print "Please select calibration file:"
      for i,tf in enumerate(files):
        print "%d  :  %s" %(i+1,tf)
      path_calib = files[int(raw_input())-1]
    else:
      path_calib = files[-1]

  path_calib = os.path.join(alignmentdir,path_calib)

  run=0
  calib = calp.CalibPars(path_calib,run)
  coord = CSPADPixCoords(calib)
  X,Y = coord.get_cspad_pix_coordinate_arrays_um  (config=None)
  x = np.concatenate(X,0)
  y = np.concatenate(Y,0)

  #calp.calibpars.setCalibParsForPath (run=0,path=path_calib)
  #cpe.cpeval.evaluateCSPadPixCoordinates (rotation, mirror)

  #ccp.cspadconfig.setCSPadConfiguration(fname, dsname, event=0)
  #quadNumsInEvent  = np.arange(4)
  #indPairsInQuads  = np.arange(32).reshape(4,-1)
  #nquads           = 4
  #nsects           = 8
#
  #nsects_in_data = max(indPairsInQuads.flatten()) + 1
  #x = np.zeros((nsects_in_data,185,388), dtype=np.float32)
  #y = np.zeros((nsects_in_data,185,388), dtype=np.float32)
 # 
  #for iq in range(len(quadNumsInEvent)) :
      #quad = int(quadNumsInEvent[iq]) # uint8 -> int
      #for segm in range(8): # loop over ind = 0,1,2,...,7
	  #ind_segm_in_arr = indPairsInQuads[quad][segm]
	  #if ind_segm_in_arr == -1 : continue
#
	  #x[ind_segm_in_arr][:] = cpe.cpeval.pix_global_x[quad][segm][:]
	  #y[ind_segm_in_arr][:] = cpe.cpeval.pix_global_y[quad][segm][:]
  return x,y

def getCsPadPixCoordinates_old(path_calib=alignmentdir+'/calib-xpp-2013-01-29', 
                           rotation=0, 
                           mirror=False):
  if not path_calib:
    allfiles = os.listdir(alignmentdir)
    files = []
    for tf in allfiles:
      if 'calib-' in tf:
        files.append(tf)
    print "Please select calibration file:"
    for i,tf in enumerate(files):
      print "%d  :  %s" %(i+1,tf)
    path_calib = files[int(raw_input())-1]
  path_calib = os.path.join(alignmentdir,path_calib)


  calp.calibpars.setCalibParsForPath (run=0, path=path_calib)
  cpe.cpeval.evaluateCSPadPixCoordinates (rotation)
  cpe.cpeval.evaluateCSPadPixCoordinatesShapedAsData()
  x,y = cpe.cpeval.getCSPadPixCoordinatesShapedAsData_um()

  #calp.calibpars.setCalibParsForPath (run=0,path=path_calib)
  #cpe.cpeval.evaluateCSPadPixCoordinates (rotation, mirror)

  #ccp.cspadconfig.setCSPadConfiguration(fname, dsname, event=0)
  #quadNumsInEvent  = np.arange(4)
  #indPairsInQuads  = np.arange(32).reshape(4,-1)
  #nquads           = 4
  #nsects           = 8
#
  #nsects_in_data = max(indPairsInQuads.flatten()) + 1
  #x = np.zeros((nsects_in_data,185,388), dtype=np.float32)
  #y = np.zeros((nsects_in_data,185,388), dtype=np.float32)
 # 
  #for iq in range(len(quadNumsInEvent)) :
      #quad = int(quadNumsInEvent[iq]) # uint8 -> int
      #for segm in range(8): # loop over ind = 0,1,2,...,7
	  #ind_segm_in_arr = indPairsInQuads[quad][segm]
	  #if ind_segm_in_arr == -1 : continue
#
	  #x[ind_segm_in_arr][:] = cpe.cpeval.pix_global_x[quad][segm][:]
	  #y[ind_segm_in_arr][:] = cpe.cpeval.pix_global_y[quad][segm][:]
  return x,y

class loadingtest(object):
  def loadcoo(self):
    return getCsPadPixCoordinates(path_calib='')

class CspadPattern(object):
  def __init__(self,Nx=1000,Ny=1000,path_calib='newest'):
    self._path_calib = path_calib
    self._path = os.path.abspath(__file__)
    self._xpx = [] 
    self._ypx = []
    self.load_coordinates()
    xmn = np.min(self.xpx)
    xmx = np.max(self.xpx)
    self.xVec = np.linspace(xmn,xmx,Nx)
    self._binxVec = np.linspace(xmn-(xmx-xmn)/Nx/2,
                                xmx+(xmx-xmn)/Nx/2,Nx+1)
    ymn = np.min(self.ypx)
    ymx = np.max(self.ypx)
    self.yVec = np.linspace(ymn,ymx,Ny)
    self._binyVec = np.linspace(ymn-(ymx-ymn)/Ny/2,
                                ymx+(ymx-ymn)/Ny/2,Ny+1)
    self.shp = np.shape(self.xpx)
    self.shpPattern = (Ny,Nx)
    xind = np.digitize(self.xpx.ravel(),self._binxVec)
    yind = np.digitize(self.ypx.ravel(),self._binyVec)
    self.binning = np.ravel_multi_index((yind-1,xind-1),(Ny,Nx))
    self.numperbin = np.bincount(self.binning)
    self.pixssz = np.array([110.,110])*1e-6
    self.bin = wrapFunc(self._bin,isPerEvt=True)


  def imageShow(self,I):
    img = self._bin(I)
    tools.imagesc(self.xVec,self.yVec,img)

  def _bin(self,I):
    Ilong = np.bincount(self.binning, weights=np.asfarray(I.ravel()))
    Ilong = Ilong/self.numperbin
    P = np.zeros(np.prod(self.shpPattern))
    P[:len(Ilong)] = Ilong
    P = np.reshape(P,self.shpPattern)
    return P

  def ind2ori(self,ind):
    if np.size(ind) == np.size(self.shpPattern):
      pass
    indbool = ind
    indbool = indbool.ravel()
    indout = indbool[np.ix_(self.binning)]
    indout = np.reshape(indout,self.shp)
    return indout

  def getRawdataIndex(self,coo):
    coo = np.asarray(coo)
    xdist = np.abs(self.xpx.ravel() - coo[0] - np.min(self.xpx.ravel()))
    ydist = np.abs(self.ypx.ravel() - coo[1] - np.min(self.ypx.ravel()))
    dist = np.sqrt(xdist**2+ydist**2)
    ind = np.argmin(dist)
    return np.unravel_index(ind,self.shp)

  def polygonmask(self,polygon):
    return polygonmask(polygon,self.xpx,self.ypx)

  def load_coordinates(self,rotation=0,mirror=0):
    self._xpx,self._ypx = getCsPadPixCoordinates(rotation=rotation,mirror=mirror,path_calib = self._path_calib)

  def load_coordinates_test(self,rotation=0,mirror=0):
      return getCsPadPixCoordinatesTT(rotation=rotation,mirror=mirror,path_calib = self._path_calib)

  def _get_xpx(self):
    if len(self._xpx)==0:
      #self._xpx = np.load(os.path.dirname(self._path) + '/' + 'cspad_x.npy')
      self._xpx,self._ypx = getCsPadPixCoordinates(rotation=0,path_calib = self._path_calib)
    return self._xpx
  xpx = property(_get_xpx)

  def _get_ypx(self):
    if len(self._ypx)==0:
      #self._ypx = np.load(os.path.dirname(self._path) + '/' + 'cspad_y.npy')
      self._xpx,self._ypx = getCsPadPixCoordinates(rotation=0,path_calib = self._path_calib)
    return self._ypx
  ypx = property(_get_ypx)

  def createWireMask(self,data=None):
    if data==None:
      import ixppy
      d = ixppy.dataset('/reg/g/xpp/data/example_data/cspad/liquid_scattering/hdf5/xpp43512-r0004.h5',['cspad'])
      data = d.cspad.rdStepData(0,range(10))
    i = np.mean(data,axis=0)
    ib = self.bin(i)
    tools.imagesc(self.xVec,self.yVec,ib)
    tools.clim_std()
    pl.draw()
    print "Roughly select wire end pairs, press middle mouse button when finished."
    wires = []
    tpos = np.array(pl.ginput(2))
    pointno = 0
    while not tpos==[]:
      tpos = np.array(tpos)
      pl.plot(tpos[:,0],tpos[:,1],'go-')
      for pno in range(2):pl.text(tpos[pno,0],tpos[pno,1],str(pointno),color='r')
      wires.append(tpos)
      pointno += 1
      tpos = pl.ginput(2)
    refwires = []
    for pos in wires:
      refwires.append(self._refinePoints(pos,ib))
    wirewidths = []
    print refwires
    for refwire in refwires:
      trad = self._getOneWire(refwire,i,rmax=1000)
      wirewidths.append(trad)
    
    masked = np.zeros(np.shape(i),dtype=bool)
    ax = self.xpx
    ay = self.ypx
    for refwire,wirewidth in zip(refwires,wirewidths):
      points = refwire
      m = np.diff(points[:,1])/np.diff(points[:,0])
      c = points[0,1] - points[0,0]*m
      m_ = -1/m
      dist_perp = (ay-m*ax-c)/np.sqrt(m**2+1)
      dist_par  = (ay-m_*ax)/np.sqrt(m_**2+1)
      masked[np.abs(dist_perp)<wirewidth] = True
    np.save('cspad_last_pixel_mask.npy',masked)


  def _refinePoints(self,points,i,refinerad=5000):
    #refine points
    newpoints = []
    reffig = pl.figure()
    for point in points:
      xind = (point[0]-refinerad <= self.xVec) & (self.xVec <= point[0]+refinerad)
      yind = (point[1]-refinerad <= self.yVec) & (self.yVec <= point[1]+refinerad)
      pl.clf()
      tools.imagesc(self.xVec[xind],self.yVec[yind],
        i[np.min(yind.nonzero()):np.max(yind.nonzero())+1,np.min(xind.nonzero()):np.max(xind.nonzero())+1])
      pl.plot(point[0],point[1],'go')
      tools.clim_std()
      pl.draw()
      newpoints.append(pl.ginput(1)[0])
    points = np.array(newpoints)
    return points

  def _getOneWire(self,points,i,rmax=1000):
    inp = 0
    while not inp=='q':
      print inp
      rmin = float(inp)
      tools.nfigure('wire radius selector')
      pl.clf()
      self._plotStripe(points,i,rmin,rmax)
      pl.draw()
      inp = raw_input('Enter radius change, q to finish: ')
      if inp is not 'q': width=int(inp)
    return width



  def _plotStripe(self,points,i,rmin,rmax):
    ax = self.xpx
    ay = self.ypx
    #find closest points to find limits
    m = np.diff(points[:,1])/np.diff(points[:,0])
    c = points[0,1] - points[0,0]*m

    m_ = -1/m
    dist_perp = (ay-m*ax-c)/np.sqrt(m**2+1)
    dist_par  = (ay-m_*ax)/np.sqrt(m_**2+1)
    distparrange = dist_par[np.abs(dist_perp)<=rmax]
    distparrange = [np.min(distparrange),np.max(distparrange)]
    i[np.abs(dist_perp)<rmin] = np.nan


    tperpvec = np.arange(-rmax,rmax,110.*np.sqrt(2)) 
    tparvec  = np.arange(distparrange[0],distparrange[1],110.*np.sqrt(2))
    perpind = np.digitize(dist_perp.ravel(),tperpvec)
    parind  = np.digitize(dist_par.ravel(),tparvec)
    binning = np.ravel_multi_index((perpind,parind),(len(tperpvec)+1,len(tparvec)+1))
    numperbin = np.bincount(binning)
    wireROI = np.bincount(binning, weights=np.asfarray(i.ravel()))
    wireROI = wireROI/numperbin
    P = np.zeros((len(tperpvec)+1) * (len(tparvec)+1))
    P[:len(wireROI)] = wireROI
    P = np.reshape(P,(len(tperpvec)+1,len(tparvec)+1))
    pl.clf()
    tools.imagesc(tparvec[:-1]+np.mean(np.diff(tparvec)),
                  tperpvec[:-1]+np.mean(np.diff(tperpvec)),
                  P[1:-1,1:-1])
    tools.clim_std()
    pl.axis('normal')
    pl.draw()

  def saveMaskDAQ(self,filename,mask=None):
    if mask==None:
      mask = self.mask
    f = file(filename,'w')
    for seg in mask:
      for line in seg:
        for el in line:
          f.write(' '+str(int(el)))
        f.write('\n')
  
  def showTileNumbers(self):
    x = [np.mean(tx) for tx in self.xpx.ravel()]
    y = [np.mean(ty) for ty in self.ypx.ravel()]
    for n,(tx,ty) in enumerate(zip(x,y)):
      plt.text(tx,ty,str(n))


  def refineCen(self,i,segments=16,Nrbins='auto',cycles=5):
    #def refineCen(x,y,i,cen,R0ring,Rsearchlims,segments=16,Nrbins='auto',cycles=5):
    self.imageShow(i)
    x = self.xpx.ravel()
    y = self.ypx.ravel()
    print "Select ring center..."
    cen = np.asarray(plt.ginput(1))[0]
    print "Select point on ring ..."
    por = np.asarray(plt.ginput(1))[0]
    R0ring = np.sqrt(np.sum((por-cen)**2))
    print "Select search limits ..."
    rslc = np.asarray(plt.ginput(2))
    Rslc = np.sqrt(np.sum((rslc-cen)**2,axis=1)).ravel()
    Rslc.sort()
    Rsearchlims = Rslc-R0ring

    azedges = np.arange(0,segments+1)*(2*np.pi/segments)-np.pi
    if Nrbins=='auto':
      rbinsz = np.max(np.diff(np.sort(x.ravel())))
      Rsearchlims = np.asarray(Rsearchlims)
    mask = None
    for cycle in range(cycles):
      redges = np.arange(R0ring+Rsearchlims[0],R0ring+Rsearchlims[1],rbinsz)
      R = np.sqrt((x-cen[0])**2+(y-cen[1])**2)
      A = np.arctan2(y-cen[1],x-cen[0])
      rd = np.digitize(R.ravel(),redges)
      azd = np.digitize(A.ravel(),azedges)
      plt.figure(5)
      plt.hold(0)
      self.imageShow(azd.reshape(x.shape))


      idx = np.ravel_multi_index(np.vstack([rd,azd]), (len(redges)+1,len(azedges)+1))
      mnl = (len(redges)+1)*(len(azedges)+1)
      res = (np.bincount(idx,weights=i.ravel(),minlength=mnl)/np.bincount(idx,minlength=mnl)).reshape(len(redges)+1,len(azedges)+1)[1:-1,1:-1]
      rvec = tools.histVecCenter(redges)
      avec = tools.histVecCenter(azedges)
      #fts = [np.polyfit(rvec[~np.isnan(tprof)],tprof[~np.isnan(tprof)],2) for tprof in res.T if np.sum(~np.isnan(tprof))>1 else np.array([np.nan]*3)]
      if False:
        fts = [np.polyfit(rvec[~np.isnan(tprof)],tprof[~np.isnan(tprof)],2) if np.sum(~np.isnan(tprof))>1 else np.array([np.nan]*3) for tprof in res.T ]
        rts = np.asarray([np.unique(np.roots(np.polyder(tp)))[0] if not np.isnan(tp).all() else np.nan for tp in fts])
      else:
	rts = []
	#rts = [tools.peakAna(rvec[~np.isnan(tprof)],tprof[~np.isnan(tprof)],3)[0] if np.sum(~np.isnan(tprof))>1 else np.nan for tprof in res.T ]
	for tprof in res.T:
	  try:
	    rts.append(np.float(tools.peakAna(rvec[~np.isnan(tprof)],tprof[~np.isnan(tprof)],3)[0]))
	  except:
	    rts.append(np.nan)
	 
	  plt.figure(100)
	  plt.hold(0)
	  plt.plot(rvec[~np.isnan(tprof)],tprof[~np.isnan(tprof)],'.-')
	  plt.axvline(rts[-1])
	  plt.draw()
	  plt.waitforbuttonpress()

      rts = np.asarray(rts)	 

      plt.figure(5)
      tools.imagesc(avec,rvec,res)
      plt.hold(1)
      if mask is None:
        plt.plot(avec,rts,'ow')
	mask = np.isnan(rts)
	mask[rts<np.min(rvec)] = True
	mask[rts>np.max(rvec)] = True
        mask[tools.maskPoints(avec,rts)] = True
      else:
        plt.plot(avec[~mask],rts[~mask],'ow')

      rts[mask] = np.nan
      plt.hold(0)

      #plt.figure(6)
      #for tp,pos,prf in zip(fts,rts,res.T):
	#plt.hold(0)
       # 
	#plt.plot(rvec,prf,'k')
	#plt.hold(1)
	#plt.axvline(pos)
	#plt.plot(rvec,np.polyval(tp,rvec),'r')
	#sleep(.1)
	#plt.draw()
      xf,yf = tools.pol2cart(avec[~np.isnan(rts)],np.squeeze(rts[~np.isnan(rts)]))
      xc,yc,R0ring,chisq = tools.fitCircle(xf+cen[0],yf+cen[1],w=np.sum(~np.isnan(res),0)[~np.isnan(rts)])
      plt.figure(1)
      plt.hold(1)
      plt.plot(xc,yc,'r+')
      plt.plot(xf+cen[0],yf+cen[1],'ow')
      plt.waitforbuttonpress()
      cen = [xc,yc]
    return cen,R0ring

  #def digitizeRadialQspace(i, 
                      #polarization=None, # no correction if None, standard XPP setup correction if True, if float polarization degree.
                      #Ephot=None, # in keV
                      #polarization=None, # no correction if None, 1 for horizontal, 0 for vertical.
                      #beamCenter=None # beam center, standard center if None
                      #detDist=None, # in m
                      #detTilt=None, # detector tilt, rotation matrix
                      #binning='pixel', # binning of data: approximate pixel resolution, linear spaced in Q (str: pixel); number of bins (int); direct entry as array of floats (in rec Angstrom).
                      #gainmap=None,             # no gain correction if None, same size as matrix otherwise, dictionary of monitor (str), array of intensity values, 
                      #liquidSheetAbsorption=None, # correction for liquid sheet absorption, keyword for different methods
                      #detectorIncidenceAngle=False, # correct for detector material transmissioni
                     #):
    ## many steps
    #icorr = i
    
    ## polarization
    ## angles: theta and phi (in det space rotating ccw, from x axis )
    #pixpos = vstack([self._xpx.ravel(),self._ypx.ravel(),zeros(np.prod(np.shape(self._xpx)))])
    #if not detTilt==None:
      #pixpos = detTilt*pixpos
    #pixpos = (pixpos.T + detDist*array([0,0,1])).T

    ## 2theta angle
    #tthe = np.arccos(dot(array([0,0,1]), pixpos) / np.apply_along_axis(np.linalg.norm,0,pixpos))
    
    ## angle to horizontal pol direction
    #polAngHor = np.arccos(dot(array([1,0,0]), pixpos) / np.apply_along_axis(np.linalg.norm,0,pixpos))
    ## angle to vertical pol direction
    #polAngVer = np.arccos(dot(array([0,1,0]), pixpos) / np.apply_along_axis(np.linalg.norm,0,pixpos))

    #polfac = 1/(polarization*cos(polAngHor)**2) + 1/((1-polarization)*cos(polAngVer)**2)
    
    ## WAXS style diffraction angles
    #azi = arcsin(pixpos[0,:]/pixpos[2,:])
    #ele = arcsin(pixpos[1,:]/np.sqrt(pixpos[0,:]**2 + pixpos[2,:]**2))



    #sint
    #sinp
    #cosp
    #polfac = 1/ (polarization*(1-(sinp*sint)**2) + (1-polarization)*(1-(cosp*sint)**2))

    ## solid angle correction
    ## delta: angle from detector surface normal
    #1/cos(theta)**3

    ## absorption correction
    ## for strainght on sheet
    #liqsheetfac = (musheet*dsheet - musheet*dsheet/cos(2*theta)) /
                  #(exp(musheet*dsheet) * -exp(-musheet*dsheet) + exp(-musheet*dsheet)*exp(cos(2*theta))
    ## solution for tilted sheet
    
    #Xsheet = (1-p) * dsheet / cos(alpha+gamma)
    #Lsheet = sqrt(Xsheet**2 + (Xsheet*tan(delta))**2)
    #Transmission = exp(-musheet*Lsheet)

    ## detector transmission length

    #Ddetfac = (1 - exp(-mudet*ddet))/
                   #(1 - exp(-mudet*ddet / cos(theta)))












    ## output
    #return binvec,Qvec
    
  #def binRadialQspace(i, 
                      #binvec,
                      #dark=None,
                      #correction=None, # correction in form of array of same size as data that is multiplied on the data before binning.
                     #):
    #i = i-dark
    ## many steps
    #icorr = i



    ## output
    #return iBinVec

#pattern = CspadPattern()
  

def noiseMap(Istack):
  return np.std(np.asfarray(Istack)/sum(Istack),axis=0)

def getNoiseMap(Istack,lims_perc=None,lims=None):
  noise = noiseMap(Istack)
  #np.shape(noise)
  if lims_perc is not None:
    lims = np.percentile(noise,lims_perc)
    tools.nfigure('Selected noise limits')
    pl.clf()
    tools.histSmart(noise.ravel()[~np.isnan(noise.ravel())],fac=200)
    pp = plt.axhspan(*lims)
    plt.gca().add_patch(pp)
    plt.draw()

  if lims==None:
    tools.nfigure('Find noise limits')
    pl.clf()
    tools.histSmart(noise.ravel()[~np.isnan(noise.ravel())],fac=200)
    #pl.gca().set_xscale('log')
    pl.draw()
    print "Select noise limits"
    lims = tools.getSpanCoordinates()
  return ~tools.filtvec(noise,lims),noise

def getDarkNoise(data,maximgs=10,xoff=None,save='',lims_perc=None,lims=None):
  if xoff is not None:
    dark_raw = xoff*data
  else:
    dark_raw = data
  lens = dark_raw.lens()
  darks = []
  remaining = maximgs
  widgets = ['Reading %d dark images: '%maximgs, pb.Percentage(), ' ', pb.Bar(),' ', pb.ETA(),'  ']
  pbar = pb.ProgressBar(widgets=widgets, maxval=maximgs).start()
  for stepNo,tlen in enumerate(lens):
    darks.append(dark_raw[stepNo,:min(tlen,remaining)][0])
    remaining -= tlen
    if remaining<1: break
    pbar.update(maximgs-remaining)
  darks = np.concatenate(darks,axis=0)
  pbar.finish()
  noisemask,noise = getNoiseMap(darks,lims_perc=lims_perc,lims=lims)
  dark = np.mean(darks,axis=0)
  if save is not None:
    def saveFile(fina,dat):
      overwrite = True
      if os.path.exists(fina):
	ip = raw_input('File %s exists, wanna overwrite? (y/n)'%fina)
	overwrite = ip=='y'
      if overwrite:
	np.save(fina,dat)
    saveFile('noisemask'+save+'.npy',noisemask)
    saveFile('noise'+save+'.npy',noise)
    saveFile('dark'+save+'.npy',dark)
  return dark,noisemask


def corrLongPix(I,fillvalues=True,BGcorrect=None):
  if BGcorrect:
    if BGcorrect<0:
      stripe = I[-BGcorrect:,:]
    elif BGcorrect>0:
      stripe = I[:BGcorrect,:]
    stripe = np.ma.masked_array(stripe,np.isnan(stripe))
    I = I-np.mean(stripe,axis=0)
  Io = np.ones((185,388+3))*np.nan
  Io[:,:388/2-1] =  I[:,:388/2-1]
  Io[:,-388/2+1:] =  I[:,388/2+1:]
  if fillvalues:
    mpp1 = I[:,388/2-1].copy()/2.5
    mpp2 = I[:,388/2].copy()/2.5
    for n in range(2):
      Io[:,388/2-1+n]  =  mpp1.copy()
      Io[:,-388/2-n] =  mpp2.copy()
    Io[:,388/2+1] = np.mean(np.vstack([mpp1,mpp2]),axis=0).transpose()
  return Io




#class cspad(object):
  #def getCommonModeFromHist(im,searchoffset=200,COMrad=3):                        
    #bins = np.arange(1000) 
    #hist,dum = histogram(im.ravel(),bins)
    #aboveoffset = (hist>searchoffset)
    #iao = list(aboveoffset).index(True)
    #imx = iao + list(diff(hist[aboveoffset])<0).index(True)
    #CM = (sum(hist[imx-COMrad:imx+COMrad+1]*bins[imx-COMrad:imx+COMrad+1])/ sum(hist[imx-COMrad:imx+COMrad+1])) 
    #return CM

  #def rdCSPAD_metrology_coordinates(fina='CXI1-Metrology-Feldkamp.xls'):
    #wb = xlrd.open_workbook(fina)
    #Pos = []
    #for shNO in range(4):
      #meanpos = []
      #tilts = []
      #isportrait = []
      #sh = wb.sheet_by_index(shNO)
      #for segNO in range(8):
        #sco = []
        #for corNO in range(4):
                #rv = sh.row_values(segNO*4+corNO+1)[1:4]
                #sco.append(rv)
        #sco = pl.array(sco)
        #sides = pl.array([sco[[0,1,2,3],0]-sco[[1,2,3,0],0],sco[[0,1,2,3],1]-sco[[1,2,3,0],1]])
        #tisportrait = diff(sum(abs(sides),axis=1))[0]>0
        #tilt,dum = cart2pol(sides[0,:],sides[1,:])
        #tilt = remainder(tilt+10,90)-10
        #meanpos.append(list(mean(sco,axis=0)/1e3))
        #tilts.append(mean(tilt))
        #isportrait.append(tisportrait)
      #Pos.append([pl.array(meanpos),tilts,isportrait])
    #return Pos


  #def CSPADassignsimple(data=0,cspadconfig=0):
    #segments = np.arange(32)
    #Pos = rdCSPAD_metrology_coordinates()
    #allxmin = []
    #allymin = []
    #allxmax = []
    #allymax = []
    #shortoffs = pl.ones(8)*185*.11/2
    #shortoffs = pl.ones(8)*185*.11/2
    #longoffs = pl.ones(8)*(2*194+3)*.11/2

    #for quad in Pos:
      #allxmin.extend(list(pl.array(quad[0])[:,0] - quad[2]*shortoffs - ~pl.array(quad[2])*longoffs))
      #allymin.extend(list(pl.array(quad[0])[:,1] - ~pl.array(quad[2])*shortoffs - quad[2]*longoffs))
      #allxmax.extend(list(pl.array(quad[0])[:,0] + quad[2]*shortoffs + ~pl.array(quad[2])*longoffs))
      #allymax.extend(list(pl.array(quad[0])[:,1] + ~pl.array(quad[2])*shortoffs + quad[2]*longoffs))
    
    #xpxtot = np.round((max(allxmax)-min(allxmin))/.11)
    #ypxtot = np.round((max(allymax)-min(allymin))/.11)
    #xtotpatt = zeros([ypxtot+1,xpxtot+1])
    #ytotpatt = xtotpatt
    #shortsegmat,longsegmat = meshgrid(np.arange(185),np.arange(194*2+3))
    #shortsegmat = shortsegmat-mean(shortsegmat)
    #longsegmat = longsegmat-mean(longsegmat)

    #for quad in Pos:
      #for segNO in np.arange(pl.shape(quad[0])[0]):
        #if quad[2][segNO]:
                #tx = quad[0][segNO,0]-min(allxmin) + shortsegmat*0.11
                #ty = quad[0][segNO,1]-min(allymin) + longsegmat*0.11
        #else:
                #tx = quad[0][segNO,0]-min(allxmin) + longsegmat.transpose()*0.11
                #ty = quad[0][segNO,1]-min(allymin) + shortsegmat.transpose()*0.11
                
        #xtotpatt[pl.ix_(np.round(ty[:,1]/.11).astype('int').tolist(),np.round(tx[1,:]/0.11).astype('int').tolist())] = tx
        #ytotpatt[pl.ix_(np.round(ty[:,1]/.11).astype('int').tolist(),np.round(tx[1,:]/0.11).astype('int').tolist())] = ty

    #return xtotpatt


def create_pixel_histogram(cspad_det,dark=None,Nmax=None):
  cc = 0
  npatt = len(cspad_dat.time[cc])
  if Nmax==None or Nmax>npatt:
    Nmax = npatt
  chunks = cspad_dat.chunks(Nmax=Nmax)

  histbins = []
  histograms = []

  for ch in chunks[cc]:
    data = ch.data
    if not dark==None:
      data = data-dark
    if histbins==[]:
      datshp = np.shape(data[0])
      Nel = np.prod(datshp)
      histbins = np.empty(datshp,dtype=np.object_)
      data = data.reshape(datshp[0],-1)
      ind = 0
      for n in range(Nel):
        ind = np.ravel_index(n,datshp)
        histbins[ind] = histEdges(data[ind])
    
    for n in range(Nel):
      ind = np.ravel_index(n,datshp)
      if histograms ==[]:
        histograms = np.empty(datshp,dtype=np.object_)
        histpgrams[ind] = np.histogram(data[ind],histbins[ind])
      else:
        histograms[ind] += np.histogram(data[ind],histbins[ind])

  return histograms,histbins


def maskEdge(shape=(185,388),offset=1,maskmid=True):
  msk = np.zeros(shape,dtype='bool')
  msk[:offset,:] = True
  msk[:,:offset] = True
  msk[-offset:,:] = True
  msk[:,-offset:] = True
  if maskmid:
    msk[:,np.floor((max(shape)-1)/2.)] = True
    msk[:,np.ceil((max(shape)-1)/2.)] = True
  return msk
    
def maskEdges(i,offset=1,maskmid=True):
  if type(i) is int:
    shp = (i,185,388)
  else:
    shp = np.shape(i)
  shpTile = shp[-2:]
  shpdet  = shp[:-2]
  mskTile = maskEdge(shape=shpTile,offset=offset,maskmid=maskmid)
  tmsk = mskTile
  iter_shpdet = list(shpdet)
  iter_shpdet.reverse()
  for N in iter_shpdet:
    tmsk = [tmsk for n in xrange(N)]

  return np.asarray(tmsk)


def createMaskUnbonded(Ntiles=32, maskDirectNeighbors=False):
  if Ntiles not in g_maskUnbonded:
    p = np.zeros([185,388],dtype=bool)
    for i in range(0,185,10):
      p[i,i] = True
      p[i,i+194] = True
      if maskDirectNeighbors:
        for dx in [-1,1]:
          p[i,i+dx] = True
          p[i,i+194+dx] = True
        for dy in [-1,1]:
          p[i+dy,i] = True
          p[i+dy,i+194] = True
    p = p.reshape(1,185,388)
    p = np.concatenate([p]*Ntiles,0)
    globals()["g_maskUnbonded"][Ntiles]  = p
  return g_maskUnbonded[Ntiles]

def getUnbondedPixelOffset(img):
  nTiles = img.shape[0]
  msk    = createMaskUnbonded(nTiles,maskDirectNeighbors=False)
  dark   = img[msk].reshape( (-1,nTiles ) ).mean(axis=0)
  return dark

def correctImageUsingUnbondedPixels(img,makeCopy=False):
  dark = getUnbondedPixelOffset(img)
  # broadcast
  if img.ndim == 3:
    dark = dark[:,np.newaxis,np.newaxis]
  elif img.ndim == 4:
    dark = dark[:,np.newaxis,np.newaxis,np.newaxis]
  if makeCopy:
    img = copy.copy(img) - dark.astype( img.dtype )
  else:
    img -= dark.astype( img.dtype )
  return img

def correctImagesUsingUnbondedPixels(imgs,makeCopy=False):
  if makeCopy:
    out = copy.copy(imgs)
  else:
    out = imgs
  for i,img in enumerate(out):
    correctImageUsingUnbondedPixels(img,makeCopy=False)
  return out


def corrAreadetNonlin_getComponents(areadet,I0,digibins=None):
  print "finding intensity intervals for component determination"
  areadet['cnl_I0'] = I0.digitize(digibins)
  areadet['cnl_I0sorted'] = I0dig.ones()*areadet.data/I0dig
  areadet.cnl_I0sorted.evaluate[:,:100]



def unbAnalysis(datastack,bins = None):
  
  msk = maskEdge()
  unb = createMaskUnbonded(1)[0]
  out = []
  for img in datastack:
    iout = []
    for tile in img:
      tunb = np.median(tile[unb])
      tile = tile[~msk]
      tint = np.median(tile)
      iout.append(np.asarray([tunb,tint]))
    out.append(np.asarray(iout))
  return np.asarray(out)

def getCommonModeFromHist(im,gainAv=30,searchRadiusFrac=.4,debug=False):                       
  im = im.ravel()
  bins = np.arange(-2*gainAv,3*gainAv)
  hst,dum = np.histogram(im.ravel(),bins)
  bins = bins[:-1]+.5
  rad = np.round(gainAv*searchRadiusFrac)
  idx = filtvec(bins,[-rad,rad])
  pk = bins[idx][hst[idx].argmax()]
  idx = filtvec(bins,[pk-rad,pk+rad])
  bins = bins[idx]
  hst = hst[idx]
  bg = np.sum(bins*hst)/np.sum(hst)
  if debug:
    nfigure('debug common mode hist correction')
    plt.clf()
    plt.plot(bins,hst)
    plt.waitforbuttonpress()
  return bg



def commonModeCorrectTile(tile,mask=None,gainAv=30,nbSwitchFactor=3,unbPx=None):
  if unbPx is None:
    unb = createMaskUnbonded(1)[0]
  else:
    unb = unbPx
  unbV = np.median(tile[unb])
  if gainAv is not None:
    if mask is None:
      mask = np.ones_like(tile,dtype=bool)
    tdat = tile[mask]
    med = np.median(tdat)
    lowexp = (med-unbV) < (nbSwitchFactor*gainAv)

    if lowexp:
      bg = getCommonModeFromHist(tdat,gainAv=gainAv)
    else:
      bg = unbV
  else:
    bg = unbV
    
  tile -= bg
  return tile,bg


def commonModeCorrectImg(img,mask=None,gainAv=30,nbSwitchFactor=3,unbPx=None):
  for tNo,tile in enumerate(img):
    if mask is None:
      tmask=None
    else:
      tmask = mask[tNo]
    tile,bg = commonModeCorrectTile(tile,mask=tmask,gainAv=gainAv,nbSwitchFactor=nbSwitchFactor,unbPx=unbPx)
  return img

commonModeCorrect = wrapFunc(commonModeCorrectImg,isPerEvt=True)

def _nanify(i,mask):
  if np.shape(mask)[-1] is 1:
    i[np.squeeze(mask,axis=-1),:] = np.nan
  else:
    i[mask,:] = np.nan
  return i

nanify = wrapFunc(_nanify,isPerEvt=False,transposeStack=True)

def correct(data,dark=None,mask=None,NpxEdgeMask=1,correctCommonMode=True,gainAv=30,nbSwitchFactor=3,Ntiles=32):
  if dark is not None:
    if np.shape(dark)[-1] is not 1: dark=dark[...,np.newaxis]
    darkcorrect = data.astype(np.float)-dark
  else:
    darkcorrect = data.astype(np.float)
  maskub = createMaskUnbonded(Ntiles)
  maskedg = maskEdges(Ntiles,offset=NpxEdgeMask)
  maskcomb = np.logical_or(mask,maskub)
  maskcomb = np.logical_or(mask,maskedg)
  if np.shape(maskcomb)[-1] is not 1: maskcomb=maskcomb[...,np.newaxis]

  if correctCommonMode:
    corr0 = commonModeCorrect(darkcorrect,mask=maskcomb,gainAv=gainAv,
      nbSwitchFactor=nbSwitchFactor,unbPx=maskub[0])
  else:
    corr0 = darkcorrect
  corr0 = nanify(1.*corr0,maskcomb)
  return corr0
  

def histOverview(data,clearFig=True):
  Nax = len(data)
  msk = maskEdge()
  unb = createMaskUnbonded(1)[0]
  fig = nfigure('Cspad histogram overview')
  if clearFig: plt.clf()
  binvec = np.arange(np.round(np.min(data.ravel())),np.round(np.max(data.ravel())),1)
  ah = []
  for n,tdat in enumerate(data):
    #tdat = commonModeCorrectTile(tdat)[0]
    unbpx = tdat[unb]
    tdat = tdat[~msk]
    if len(ah)==0:
      ah.append(plt.subplot(8,4,n+1))
    else:
      ah.append(plt.subplot(8,4,n+1,sharex=ah[0]))
    h = np.histogram(tdat,bins=binvec)
    lh = plt.step(binvec[:-1],h[0],where='pre')
    lh = lh[0]
    plt.axvline(np.median(unbpx),color=lh.get_color())
    plt.axvline(np.median(tdat),linestyle='--',color=lh.get_color())
    plt.text(.5,.8,str(n),horizontalalignment='center',transform=ah[-1].transAxes)
    #plt.title(str(n))
  fig.subplots_adjust(hspace=0)


class corrNonLin(object):
  def __init__(self,data):
    self.data = data
    self.refDataFilter = 1

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


def genRayonixCoordinates(binningNo=2):
  """make 2 2d arrays with x and z coordinates for the rayonix."""
  edges = np.arange(0,3841,binningNo)*.44
  centers = edges[:-1]+np.diff(edges)
  return np.meshgrid(centers,centers)




