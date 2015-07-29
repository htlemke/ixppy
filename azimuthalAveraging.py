import numpy as np
import pylab as plt
import sys
from scipy import hypot,arcsin,arccos
import os
import time
import h5py
from scipy.interpolate import griddata
import utilities as util
from toolsVecAndMat import polyFit
g_az = None

def displayimg(img,**kwargs):
  plt.imshow(img.transpose(),origin="lower",**kwargs)

class azimuthalAveraging:
  def __init__(self,mask,xcen,ycen,x=None,y=None,pixelsize=100e-6,d=100e-3,tx=0,ty=0,
      qbin=5e-3,lam=1,verbose=0,Pplane=1,phibin=0.1,img=None,report_file="auto"):
    """ pixelsize = pixel size (in m)
        tx,ty = angle of detector normal with respect to incoming beam (in deg)
                zeros are for perpendicular configuration
        x,y   = pixel coordinate (if None automatic calculation)
        qbin = rebinning q 
        phibin = bin in azimuthal angle (used for polar plot
        Pplane = Polarization (1 = horizontal, 0 = vertical)
        d     = distance of center of detector to sample (in m)
        lam   = wavelength in Ang
        img is used only for displaying corrections
        TODO: Polarization correction
    """
    self.verbose=verbose
    mask = np.asarray(mask,dtype=np.bool)
    tx = np.deg2rad(tx)
    ty = np.deg2rad(ty)
    xcen = float(xcen)
    ycen = float(ycen)
    # equations based on J Chem Phys 113, 9140 (2000) [logbook D30580, pag 71]
    (A,B,C) = (-np.sin(ty)*np.cos(tx),-np.sin(tx),-np.cos(ty)*np.cos(tx))
    (a,b,c) = (xcen+d/pixelsize*np.tan(ty),float(ycen)-d/pixelsize*np.tan(tx),d/pixelsize)
    #(a,b,c) = (xcen,ycen,d/pixelsize)
    #A=-0.884
    #B=0
    #C=-0.467
    #a=997.
    #b=1161.
    #c=573.
    #print "A,B,C",A,B,C
    #print "a,b,c",a,b,c
    self.xcen = xcen
    self.ycen = ycen
    mshape = mask.shape
    if (x is None) or (y is None):
      (x,y)=np.mgrid[0:mshape[0],0:mshape[1]]
      x += (+0.5)
      y += (+0.5)
    r  = np.sqrt( (x-a)**2+(y-b)**2+c**2)
    self.msg("calculating theta...",cr=0)
    matrix_theta = np.arccos( (A*(x-a)+B*(y-b)-C*c )/r )
    self.matrix_theta = matrix_theta
    self.msg("...done")
    #(a,b,c) = (xcen,ycen,d/pixelsize)
    # r  = n.sqrt( (x-a)**2+(y-b)**2+c**2)
    self.msg("calculating phi...",cr=0)
    matrix_phi   = np.arccos( ((A**2+C**2)*(y-b)-A*B*(x-a)+B*C*c )/ \
        np.sqrt((A**2+C**2)*(r**2-(A*(x-a)+B*(y-b)-C*c)**2)))
    idx = (y>ycen) & (np.isnan(matrix_phi))
    matrix_phi[idx] = 0
    idx = (y<ycen) & (np.isnan(matrix_phi))
    matrix_phi[idx] = np.pi
    idx = (x<xcen)
    matrix_phi[idx] = (np.pi-matrix_phi[idx])+np.pi
#    matrix_phi[idx] = temp+n.pi
    self.matrix_phi = matrix_phi
    self.msg("...done")
    Pout   = 1-Pplane
    self.msg("calculating pol matrix...",cr=0)
    pol = Pout*(1-(np.sin(matrix_phi)*np.sin(matrix_theta))**2)+\
        Pplane*(1-(np.cos(matrix_phi)*np.sin(matrix_theta))**2)
    self.msg("... done")
    self.pol=pol
    self.mask=mask
    theta_max = np.nanmax(matrix_theta[mask])
    self.msg("calculating digitize")
    self.matrix_q = 4*np.pi/lam*np.sin(self.matrix_theta/2)
    q_max = np.nanmax(self.matrix_q[mask])
    qbin = np.array(qbin)
    if qbin.size==1:
      self.qbins = np.arange(0,q_max+qbin,qbin)
    else:
      self.qbins = qbin
    self.q = (self.qbins[0:-1]+self.qbins[1:])/2
    self.theta = 2*np.arcsin(self.q*lam/4/np.pi)
    self.nq = self.q.size
    self.idxs  = np.digitize(self.matrix_q.ravel(),self.qbins)
    last_idx = self.idxs.max()
    self.idxs[~mask.ravel()] = 0; # send the masked ones in the first bin
    #print "last index",last_idx
    self.msg("...done")
    self.phi  = np.arange(0,2*np.pi+phibin,phibin)+phibin/2
    # include geometrical corrections
    self.msg("calculating normalization...",cr=0)
    geom  = (d/r) ; # pixels are not perpendicular to scattered beam
    geom *= (d/r**2); # scattered radiation is proportional to 1/r^2
    self.geom = geom
    self.geom /= self.geom.max()
    self.correction = self.geom*self.pol
    #norm = n.bincount(self.idxs,self.geom.ravel(),minlength=self.nq)
    self.Npixel = np.bincount(self.idxs,minlength=self.nq)
    self.Npixel = self.Npixel[:self.nq]
    #self.correction1D  =np.bincount(self.idxs,np.ravel(1/self.correction),minlength=self.nq)
    #self.correction1D  =self.correction1D[:self.nq]/self.Npixel
    self.norm=self.Npixel;#/self.correction1D
    self.header  = "# Parameters for data reduction\n"
    self.header += "# xcen,ycen = %.2f %.2f\n" % (xcen,ycen)
    self.header += "# pixel size = %.2g m\n" % (pixelsize)
    self.header += "# sample det distance = %.4f m\n" % (d)
    self.header += "# wavelength = %.4f Ang\n" % (lam)
    self.header += "# detector angles x,y = %.3f,%.3f deg\n" % (np.rad2deg(tx),np.rad2deg(ty))
    self.header += "# fraction of inplane pol %.3f\n" % (Pplane)
    if isinstance(qbin,float):
      self.header += "# q binning : %.3f Ang-1\n" % (qbin)
    return 
    # prepare report
    if (img is None): img=np.ones_like(mask)
    plt.interactive(0)
    plt.figure(figsize=(8*2, 6*2),dpi=150)
    plt.subplot("231",title="Polarization")
    plt.imshow(self.pol)
    plt.colorbar()
    plt.subplot("232",title="Geometrical")
    plt.imshow(self.geom)
    plt.colorbar()
    plt.subplot("233",title="Geometrical+Pol")
    plt.imshow(self.correction)
    plt.colorbar()
    plt.subplot("234",title="Raw image")
    plt.imshow(img*mask)
    plt.colorbar()
    plt.subplot("235",title="Corrected image")
    plt.imshow(img/self.correction*mask)
    plt.colorbar()
#    plt.show()
    if (report_file == "auto"):
      report_file="azimuthal_averaging_info.png"
    plt.savefig(report_file)
    self.msg("...done")

  def msg(self,s,cr=True):
    if (self.verbose):
      if (cr):
        print s 
      else:
        print s,
    sys.stdout.flush()

  def displaypolar(self,img,withCorrection=True):
    self.idxsphi  = np.digitize(self.matrix_phi.ravel(),self.phi)
    if (withCorrection):
      ii = img.ravel()/self.correction.ravel()
    else:
      ii = img.ravel()
    t = self.theta
    p = self.phi
    v = np.zeros( (len(self.theta),len(self.phi) ) )
    N = np.zeros( (len(self.theta),len(self.phi) ) )
    for i in range(len(ii)):
        try:
          it = self.idxs[i]
          ip = self.idxsphi[i]
          v[it,ip] += ii[i]
          N[it,ip] += 1
        except:
          #print "skipping",i
          pass
        if (i%100000==0):
          print "done %.1f per cent" % (float(i)/len(ii)*100)
    idx = N>0
    v[idx] /= N[idx]
    #plt.imshow(v,extent=[p[0],p[-1],t[0],t[-1]],origin="bottom",interpolation="none")
    plt.subplot("221")
    plt.imshow(v)
    plt.axis('tight')
    plt.colorbar()
    plt.subplot("222")
    plt.plot(self.phi,v[300,:])
    plt.show()
    return v

  def doAzimuthalAveraging(self,img):
    #img /= self.correction
    t0=time.time()
    I=np.bincount(self.idxs,img.ravel()/self.correction.ravel(),minlength=self.nq)
    I=I[:self.nq]
    self.sig = np.sqrt(I)/self.norm
    #I2=n.bincount(self.idxs,n.ravel(img*img),minlength=self.nq)
    self.I = I/self.norm
    #self.sig = n.sqrt((I2/self.norm-self.I**2)/self.norm)
    #idx = self.Npixel>20
    #self.sig[~idx]=self.sig2[~idx]
    return self.I

  def saveChiFile(self,fname):
    header = "q(Ang-1) I sig"
    #mc.writev(fname,self.q,n.vstack((self.I,self.sig)),header=header)
    #n.savetxt(fname,n.vstack((self.q,self.I,self.sig)),header=header)
    np.savetxt(fname,np.transpose(np.vstack((self.q,self.I,self.sig))),
      fmt=["%.3f","%.4f","%.4f"])


def readMar(f):
  from fabio import marccdimage as mccd
  print "reading image",f
  i=mccd.marccdimage()
  i.read(f)
  # conversion to float32 is fastest ... (but then in bincount is converted to float anyway)
  #return i.data
  return i.data.astype(np.float)

def doFolderOrFiles(folderNameOrFileList,
    skipFirst=0,
    forceChi=False,
    psize=100e-6,
    d = 0.1,
    xcen = 1024,
    ycen = 1024,
    lam = 1.,
    qbin = 0.01,
    tx = 0,
    ty = 0,
    x=None,
    y=None,
    folderOut=None,
    mask = None,
    waitForFiles=True,
    imageReader = readMar,
    ccdfilenames="*.mccd"
    ):
  """ Perform azimuthal averaging.
    folderNameOrFileList: can be either a folder (where the *.mccd will be found)
                          or a file list
    skipFirst: to skip the first files (but careful to order...)
    forceChi : calculate chi file even if output chi is present
    d = sample detector distance 
    pixelsize = pixel size (in m)
    xcen,ycen = center of the image
    tx,ty = angle of detector normal with respect to incoming beam (in deg)
        zeros are for perpendicular configuration
    qbin = rebinning q (spacing or list)
    Pplane = Polarization (1 = horizontal, 0 = vertical)
    d     = distance of center of detector to sample (in m)
    lam   = wavelength in Ang
    folderOut : if None, same as folderNameOrFileList
    imageReader : function that takes a name and return the intensity matrix
    ccdfilenames : pattern to look for files
  """
  az = None
  t_start = util.now()
  while ((az is None) or (waitForFiles)):
    if (os.path.isdir(folderNameOrFileList)):
      f=os.popen("ls -1 %s/%s" % (folderNameOrFileList,ccdfilenames))
      temp=f.readlines()
      f.close()
      files = []
      for t in temp:
        files.append(t.strip())
      files = files[skipFirst:]
      if folderOut is None: folderOut = folderNameOrFileList
    else:
      if type(folderNameOrFileList) != list:
        folderNameOrFileList = (folderNameOrFileList,); # if we pass single file
      files=folderNameOrFileList[skipFirst:]
      if folderOut is None: folderOut = os.path.dirname(folderNameOrFileList[0])
    if (az is None):
      if not os.path.exists(folderOut): os.makedirs(folderOut)
      t0=time.time()
      f = files[0]
      iData = imageReader(f)
      i0_mask = iData != 0
      fname = folderOut+"/"+"azimuthal_averaging_info.png"
      az = azimuthal_averaging(mask&i0_mask,xcen,ycen,pixelsize=psize,x=x,y=y,d=d,
          lam=lam,qbin=qbin,img=iData,report_file=fname)
      #az.displaypolar(iData)
      print "Time needed for inizialization %.2f"%(time.time()-t0)
      fname = folderOut+"/"+"azimuthal_averaging_info.txt"
      finfo=open(fname,"w")
      finfo.write(az.header)
      finfo.close()
    if (len(files) == 0):
      print "Done %d files I could find, waiting for new files (%s)" % (skipFirst,util.now())
      time.sleep(10)
    t0=time.time()
    t_save=0.
    t_read=0.
    t_az=0.
    skip = 0
    flist = []
    data = []
    err  = []
    for f in files:
      fout = os.path.splitext(os.path.basename(f))[0]
      fout = folderOut + "/" + fout + ".chi"
      if (os.path.exists(fout) and not forceChi):
        skip += 1
        continue
      else:
        t1 = time.time()
        iData = imageReader(f)
        t_read += (time.time()-t1)
        t1 = time.time()
        az.do_azimuthal_averaging(iData)
        t_az += (time.time()-t1)
        t1 = time.time()
        az.saveChiFile(fout)
        t_save += (time.time()-t1)
        flist.append(f)
        data.append(az.I)
        err.append(az.sig)
    if ((len(files)-skip)!=0):
      nfiles = len(files)
      s="Time needed for %d files: %.2f ms/file"%(nfiles,(time.time()-t0)/nfiles*1e3)
      ttot = t_read+t_save+t_az
      s+="\n"+ "Fraction of time to read,calc,save     : %.2f,%.2f,%.2f" % (t_read/ttot,t_az/ttot,t_save/ttot)
      s+="\n"+ "Time per file    to read,calc,save (ms): %.2f,%.2f,%.2f" % \
          (t_read/nfiles*1e3,t_az/nfiles*1e3,t_save/nfiles*1e3)
      print s
      finfo=open(fname,"a")
      finfo.write(s)
      finfo.close()
      hname = folderOut + "/" + folderOut.rstrip("/").split("/")[-1]+".h5"
      t_end = util.now()
      if (~forceChi & (os.path.exists(hname))):
        hchi = h5py.File(hname,"r")
        flist = [hchi["flist"].value,flist]
        data = [hchi["data"].value,data]
        err  = [hchi["err"].value,err]
        hchi.close()
      hchi = h5py.File(hname,"w")
      hchi.attrs["time_start"]=t_start
      hchi.attrs["time_end"]=t_end
      hchi.attrs["info"] = az.header
      hchi.attrs["time_bench"] = s
      hchi.create_dataset("flist",data=flist)
      hchi.create_dataset("data",data=data)
      hchi.create_dataset("err",data=err)
      hchi.create_dataset("q",data=az.q)
      hchi.create_dataset("theta",data=az.theta)
      hchi.create_dataset("Npixel",data=az.Npixel)
      hchi.close()
    skipFirst += len(files)-skip
  return az


def _doImages(listOfImgs,azObj):
  nImg=len(listOfImgs)
  dataI=np.empty( (nImg,azObj.nq) )
  dataE=np.empty( (nImg,azObj.nq) )
  for i in range(nImg):
    img = listOfImgs[i]
    azObj.do_azimuthal_averaging(img)
    dataI[i,:] = azObj.I
    dataE[i,:] = azObj.sig
  return dataI,dataE

def doImages(listOfImgs,
    psize=100e-6,
    d = 0.1,
    xcen = 1024,
    ycen = 1024,
    lam = 1.,
    qbin = 0.01,
    tx = 0,
    ty = 0,
    x=None,
    y=None,
    folderOut="./",
    mask = None,
    force=False,
    nJobs = 4,
    hdf5out = None
    ):
  """ Perform azimuthal averaging.
    listOfImgs: list of images previously read
    d = sample detector distance 
    pixelsize = pixel size (in m)
    xcen,ycen = center of the image
    tx,ty = angle of detector normal with respect to incoming beam (in deg)
        zeros are for perpendicular configuration
    qbin = rebinning q (spacing or list)
    Pplane = Polarization (1 = horizontal, 0 = vertical)
    d     = distance of center of detector to sample (in m)
    lam   = wavelength in Ang
    force = if True force reinizialization
    hdf5out = if not None, use as outfile name
    folderOut : 
  """
  print "NJOBS",nJobs
  t_start = util.now()
  t0=time.time()
  if (len(listOfImgs) == 0):
    return
  if (g_az is None):
    if (mask is None):
      mask=np.ones_like(listOfImgs[0],dtype=np.bool)
    t0=time.time()
    fname = folderOut+"/"+"azimuthal_averaging_info.png"
    az = azimuthal_averaging(mask,xcen,ycen,pixelsize=psize,x=x,y=y,d=d,
        lam=lam,qbin=qbin,report_file=fname)
    print "Time needed for inizialization %.2f"%(time.time()-t0)
    globals()["g_az"] = az
  else:
    az = g_az
  fname = folderOut+"/"+"azimuthal_averaging_info.txt"
  t0=time.time()
  nq = az.nq
  nImg = len(listOfImgs)
  if (nJobs > 1):
    import jobManager
    ag = jobManager.myAgent(nJobs=nJobs,parallel="thread")
    #ag = jobManager.myAgent(nMax=nJobs,parallel="process")
    n = int(np.ceil(float(nImg)/nJobs))
    for i in range(nJobs):
      m=i*n
      M=(i+1)*n
      M=min(M,nImg)
      ag.addJob( _doImages,(listOfImgs[m:M],az) )
    ag.waitUntilAllDone(update=0.05)
    #time.sleep(10)
    dataI=np.vstack ( [x[0] for x in ag.data] )
    #dataI=np.reshape(dataI, (nImg,az.nq) )
    dataE=np.vstack ( [x[1] for x in ag.data] )
    #dataE=np.reshape(dataE, (nImg,az.nq) )
  else:
    dataI,dataE=_doImages(listOfImgs,az)
  t_end = util.now()
  s="Time needed for %d images: %.2f ms/img"%(nImg,(time.time()-t0)/nImg*1e3)
  print s
  finfo=open(fname,"a")
  finfo.write(s)
  finfo.close()
  if hdf5out is not None:
    hname = folderOut + "/" + hdf5out
    hchi = h5py.File(hname,"w")
    hchi.attrs["time_start"]=t_start
    hchi.attrs["time_end"]=t_end
    hchi.attrs["info"] = az.header
    hchi.attrs["time_bench"] = s
    hchi.create_dataset("data",data=dataI)
    hchi.create_dataset("err",data=dataE)
    hchi.create_dataset("q",data=az.q)
    hchi.create_dataset("theta",data=az.theta)
    hchi.create_dataset("Npixel",data=az.Npixel)
    hchi.close()
  return dataI,dataE,az

class azimuthalBinning:
  def __init__(self,x,y,xcen,ycen,d=100e-3,mask=None,gainImg=None,darkImg=None,tx=0,ty=0, qbin=5e-3,lam=1,\
        ADU_per_photon = 1.,Pplane=0,phibin=0.1,phiBins=1,img=None,verbose=0,report_file=None):
    """ 
        correctedImage = (Image-darkImg)/gainImg/geom_correction/pol_correction
        x,y      = pixel coordinate (1D array each); note: they should be the center of the pixels
        xcen,ycen = center beam position
        tx,ty = angle of detector normal with respect to incoming beam (in deg)
                zeros are for perpendicular configuration
        darkImg  = darkImage to subbract
        ADU_per_photon : used to estimate errors
        qbin = rebinning q 
        phibin = bin in azimuthal angle (used for polar plot
        Pplane = Polarization (1 = horizontal, 0 = vertical)
        d     = distance of center of detector to sample (in m)
        lam   = wavelength in Ang
        img is used only for displaying corrections
    """
    # save parameters for later use
    self.gainImg=gainImg
    self.darkImg=darkImg
    if mask is not None: mask = np.asarray(mask,dtype=np.bool)
    self.mask=mask
    self.verbose=verbose
    self.ADU_per_photon=ADU_per_photon

    tx = np.deg2rad(tx)
    ty = np.deg2rad(ty)
    xcen = float(xcen)
    ycen = float(ycen)
    # equations based on J Chem Phys 113, 9140 (2000) [logbook D30580, pag 71]
    (A,B,C) = (-np.sin(ty)*np.cos(tx),-np.sin(tx),-np.cos(ty)*np.cos(tx))
    (a,b,c) = (xcen+d*np.tan(ty),float(ycen)-d*np.tan(tx),d)
    self.xcen = xcen
    self.ycen = ycen
    mshape = x.shape

    r  = np.sqrt( (x-a)**2+(y-b)**2+c**2)
    self.r = r
    self.d = d
    
    self.msg("calculating theta...",cr=0)
    matrix_theta = np.arccos( (A*(x-a)+B*(y-b)-C*c )/r )
    self.matrix_theta = matrix_theta
    self.msg("...done")
    
    self.msg("calculating phi...",cr=0)
    matrix_phi   = np.arccos( ((A**2+C**2)*(y-b)-A*B*(x-a)+B*C*c )/ \
        np.sqrt((A**2+C**2)*(r**2-(A*(x-a)+B*(y-b)-C*c)**2)))
    idx = (y>ycen) & (np.isnan(matrix_phi))
    matrix_phi[idx] = 0
    idx = (y<ycen) & (np.isnan(matrix_phi))
    matrix_phi[idx] = np.pi
    idx = (x<xcen)
    matrix_phi[idx] = (np.pi-matrix_phi[idx])+np.pi
#    matrix_phi[idx] = temp+n.pi
    self.matrix_phi = matrix_phi
    self.msg("...done")

    self.msg("calculating pol matrix...",cr=0)
    Pout   = 1-Pplane
    pol = Pout*(1-(np.sin(matrix_phi)*np.sin(matrix_theta))**2)+\
        Pplane*(1-(np.cos(matrix_phi)*np.sin(matrix_theta))**2)

    self.msg("... done")
    self.pol=pol
    theta_max = np.nanmax(matrix_theta[~mask])
    self.msg("calculating digitize")
    self.nphi = phiBins
    #if phiBins > 1:
    phiint = 2*np.pi/phiBins
    pbm = self.matrix_phi + phiint/2
    pbm[pbm>=2*np.pi] -= 2*np.pi
    self.phiVec = np.linspace(0,2*np.pi+np.spacing(np.min(pbm)),phiBins+1)

    self.idxphi = np.digitize(pbm.ravel(),self.phiVec)-1
    self.matrix_q = 4*np.pi/lam*np.sin(self.matrix_theta/2)
    q_max = np.nanmax(self.matrix_q[~mask])
    qbin = np.array(qbin)
    if qbin.size==1:
      self.qbins = np.arange(0,q_max+qbin,qbin)
    else:
      self.qbins = qbin
    self.q = (self.qbins[0:-1]+self.qbins[1:])/2
    self.theta = 2*np.arcsin(self.q*lam/4/np.pi)
    self.nq = self.q.size
    self.idxq  = np.digitize(self.matrix_q.ravel(),self.qbins)-1
    last_idx = self.idxq.max()
    self.idxq[mask.ravel()] = 0; # send the masked ones in the first bin

    # 2D binning!
    self.Cake_idxs = np.ravel_multi_index((self.idxphi,self.idxq),(self.nphi,self.nq))
    self.Cake_idxs[mask.ravel()] = 0; # send the masked ones in the first bin
    #print "last index",last_idx
    self.msg("...done")
    #self.phi  = np.arange(0,2*np.pi+phibin,phibin)+phibin/2
    self.phi  = self.phiVec[:-1]
    # include geometrical corrections
    geom  = (d/r) ; # pixels are not perpendicular to scattered beam
    geom *= (d/r**2); # scattered radiation is proportional to 1/r^2
    self.msg("calculating normalization...",cr=0)
    self.geom = geom
    self.geom /= self.geom.max()
    self.correction = self.geom*self.pol
    self.Npixel = np.bincount(self.idxq,minlength=self.nq); self.Npixel = self.Npixel[:self.nq]
    self.norm   = self.Npixel
    self.Cake_Npixel = np.bincount(self.Cake_idxs,minlength=self.nq*self.nphi)
    #self.Cake_Npixel = self.Npixel[:self.nq*self.nphi]
    self.Cake_norm=np.reshape(self.Cake_Npixel,(self.nphi,self.nq));#/self.correction1D
    #self.correction1D  =self.correction1D[:self.nq]/self.Npixel
    self.header  = "# Parameters for data reduction\n"
    self.header += "# xcen,ycen = %.2f m %.2f m\n" % (xcen,ycen)
    self.header += "# sample det distance = %.4f m\n" % (d)
    self.header += "# wavelength = %.4f Ang\n" % (lam)
    self.header += "# detector angles x,y = %.3f,%.3f deg\n" % (np.rad2deg(tx),np.rad2deg(ty))
    self.header += "# fraction of inplane pol %.3f\n" % (Pplane)
    if isinstance(qbin,float):
      self.header += "# q binning : %.3f Ang-1\n" % (qbin)
    return 
    if report_file is None:
      return
    else:
      # prepare report
      if (img is None): img=np.ones_like(mask)
      plt.interactive(0)
      plt.figure(figsize=(8*2, 6*2),dpi=150)
      plt.subplot("231",title="Polarization")
      plt.imshow(self.pol)
      plt.colorbar()
      plt.subplot("232",title="Geometrical")
      plt.imshow(self.geom)
      plt.colorbar()
      plt.subplot("233",title="Geometrical+Pol")
      plt.imshow(self.correction)
      plt.colorbar()
      plt.subplot("234",title="Raw image")
      plt.imshow(img*mask)
      plt.colorbar()
      plt.subplot("235",title="Corrected image")
      plt.imshow(img/self.correction*mask)
      plt.colorbar()
#      plt.show()
      if (report_file == "auto"):
        report_file="azimuthal_averaging_info.png"
      plt.savefig(report_file)
    self.msg("...done")

  def msg(self,s,cr=True):
    if (self.verbose):
      if (cr):
        print s 
      else:
        print s,
    sys.stdout.flush()

  def displayCake(self,img,applyCorrection=True):
    ii =  self.doCake(img,applyCorrection=applyCorrection)
    plt.subplot("221")
    plt.imshow(ii)
    plt.axis('tight')
    plt.colorbar()
    plt.subplot("222")
    plt.plot(self.phi,ii[300,:])
    plt.show()
    return ii

  def doAzimuthalAveraging(self,img,applyCorrection=True):
    if self.darkImg is not None: img-=self.darkImg
    if self.gainImg is not None: img/=self.gainImg
    if applyCorrection:
      I=np.bincount(self.idxq, weights = img.ravel()/self.correction.ravel(), minlength=self.nq); I=I[:self.nq]
    else:
      I=np.bincount(self.idxq, weights = img.ravel()                        , minlength=self.nq); I=I[:self.nq]
    self.sig = np.sqrt(1./self.ADU_per_photon)*np.sqrt(I)/self.norm
    self.I = I/self.norm
    return self.I


  def doCake(self,img,applyCorrection=True):
    if self.darkImg is not None: img-=self.darkImg
    if self.gainImg is not None: img/=self.gainImg
    if applyCorrection:
      I=np.bincount(self.Cake_idxs, weights = img.ravel()/self.correction.ravel(), minlength=self.nq*self.nphi); I=I[:self.nq*self.nphi]
    else:
      I=np.bincount(self.Cake_idxs, weights = img.ravel()                        , minlength=self.nq*self.nphi); I=I[:self.nq*self.nphi]
    I = np.reshape(I,(self.nphi,self.nq))
    self.sig = 1./np.sqrt(self.ADU_per_photon)*np.sqrt(I)/self.Cake_norm  # ??? where comes this sqrt from? Ah I see...
    self.Icake = I/self.Cake_norm
    return self.Icake

  def saveChiFile(self,fname):
    header = "q(Ang-1) I sig"
    #mc.writev(fname,self.q,n.vstack((self.I,self.sig)),header=header)
    #n.savetxt(fname,n.vstack((self.q,self.I,self.sig)),header=header)
    np.savetxt(fname,np.transpose(np.vstack((self.q,self.I,self.sig))),
      fmt=["%.3f","%.4f","%.4f"])


def sepS0S2(D,azi):
  #Dma = np.ma.masked_array(D,np.isnan(D))
  p2 = 0.5*(3.*np.cos(azi)**2-1)
  if np.isnan(D).any():
    res = np.zeros([2,np.shape(D)[1]])
    for n,prof in enumerate(D.T):
      idx = ~np.isnan(prof)
      if sum(idx)<5:
        res[:,n] = np.array([np.nan,np.nan]) 
      else:
        res[:,n] = np.polyfit(p2[idx],prof[idx],1)
  else:
    res = polyFit(p2,D,order=1)
  res[0] *=-1
  return res


def test():
  mask=np.ones( (2000,2000) )
  az=azimuthal_averaging(mask,-80,1161,pixelsize=82e-6,d=4.7e-2,tx=0,ty=90-28.,thetabin=1e-1,lam=1,verbose=1)
  plt.subplot("121")
  displayimg(np.rad2deg(az.matrix_theta))
  print az.matrix_theta.min()
  plt.clim(0,180)
  plt.colorbar()
  plt.subplot("122")
  displayimg(np.rad2deg(az.matrix_phi))
  print az.matrix_phi.min(),az.matrix_phi.max()
  plt.colorbar()
  plt.clim(0,360)






if (__name__=="__main__"):
  test()
