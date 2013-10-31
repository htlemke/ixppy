import pylab as pl
import numpy as np
import os,sys
from ixppy import tools,wrapFunc
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

def getNoiseMap(Istack,lims=None):
  noise = noiseMap(Istack)
  #np.shape(noise)
  if lims==None:
    tools.nfigure('Find noise limits')
    pl.clf()
    tools.histSmart(noise.ravel()[~np.isnan(noise.ravel())])
    print "Select noise limits"
    lims = tools.getSpanCoordinates()
  return tools.filtvec(noise,lims)


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



def getCommonModeFromHist(im,searchoffset=200,COMrad=3):                        
  bins = np.arange(-100,100) 
  hist,dum = np.histogram(im.ravel(),bins)
  aboveoffset = (hist>searchoffset)
  iao = list(aboveoffset).index(True)
  imx = iao + list(np.diff(hist[aboveoffset])<0).index(True)
  CM = (1.*np.sum(hist[imx-COMrad:imx+COMrad+1]*bins[imx-COMrad:imx+COMrad+1])/ np.sum(hist[imx-COMrad:imx+COMrad+1])) 
  return CM

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



    
    








