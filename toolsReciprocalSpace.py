import numpy as np
from numpy import sin,cos,deg2rad,rad2deg,tan,pi
from toolsVarious import iterfy,dict2class
from toolsVecAndMat import rotmat3D
from toolsConstsAndConv import lam2E,E2lam




class crystal(object):
  def __init__(self,uc={'a':2*pi,'b':2*pi,'c':2*pi,'alpha':90,'beta':90,'gamma':90},crystalRotMat=np.asmatrix(np.eye(3))):
    """Unit cell angles in degrees"""
    if type(uc)==dict:
      self.uc = dict2class(uc)
    elif type(uc)==tuple:
      self.uc = self._get_uc_from_tuple(uc)
    else:
      self.uc = uc
    self.crystalRotMat = crystalRotMat

  def ucvectors(self):
    uc = self.uc
    av = np.array([uc.a,0,0])
    bv = np.array([uc.b*cos(deg2rad(uc.gamma)),uc.b*sin(deg2rad(uc.gamma)),0])
    cx = uc.c*cos(deg2rad(uc.beta))          
    cy = 1./sin(deg2rad(uc.gamma))*(uc.c*cos(deg2rad(uc.alpha))-cx*cos(deg2rad(uc.gamma)))
    cz = sqrt(uc.c**2-cx**2-cy**2)
    cv = np.array([cx, cy, cz]);
    R = self.crystalRotMat
    av = np.asmatrix(R)*np.asmatrix(av).ravel().transpose()
    bv = np.asmatrix(R)*np.asmatrix(bv).ravel().transpose()
    cv = np.asmatrix(R)*np.asmatrix(cv).ravel().transpose()
    return np.asarray(av.transpose()),np.asarray(bv.transpose()),np.asarray(cv.transpose())

  def ucvectors0(self):
    uc = self.uc
    av = np.array([uc.a,0,0])
    bv = np.array([uc.b*cos(deg2rad(uc.gamma)),uc.b*sin(deg2rad(uc.gamma)),0])
    cx = uc.c*cos(deg2rad(uc.beta))          
    cy = 1./sin(deg2rad(uc.gamma))*(uc.c*cos(deg2rad(uc.alpha))-cx*cos(deg2rad(uc.gamma)))
    cz = sqrt(uc.c**2-cx**2-cy**2)
    cv = np.array([cx, cy, cz]);
    R = self.crystalRotMat
    return np.asarray(av.transpose()),np.asarray(bv.transpose()),np.asarray(cv.transpose())

  def reclattvectors(self):
    av,bv,cv = self.ucvectors()
    av=av.ravel();bv=bv.ravel();cv=cv.ravel();
    Vc=dot(av, cross(bv,cv));
    avr=2*pi/Vc*cross(bv, cv);
    bvr=2*pi/Vc*cross(cv, av);
    cvr=2*pi/Vc*cross(av, bv);
    return avr,bvr,cvr

  def reclattvectors0(self):
    av,bv,cv = self.ucvectors0()
    av=av.ravel();bv=bv.ravel();cv=cv.ravel();
    Vc=dot(av, cross(bv,cv));
    avr=2*pi/Vc*cross(bv, cv);
    bvr=2*pi/Vc*cross(cv, av);
    cvr=2*pi/Vc*cross(av, bv);
    return avr,bvr,cvr
 
  def getQpQn(self,Q):
    Qnorm = self.getQnorm(Q)
    Qn = dot(Q,np.array([0,0,1]))
    Qp = sqrt(Qnorm**2-Qn**2)
    return Qp,Qn

  def getQnorm(self,Q):
    Qnorm = numpy.linalg.norm(Q)
    return Qnorm

  def getQhkl(self,hkl):
    avr,bvr,cvr = self.reclattvectors()
    Q = hkl[0]*avr+hkl[1]*bvr+hkl[2]*cvr
    return Q

  def getQhkl0(self,hkl):
    avr,bvr,cvr = self.reclattvectors0()
    Q = hkl[0]*avr+hkl[1]*bvr+hkl[2]*cvr
    return Q
  #def _transform_uc(self,uc):
    #if type(uc) is not dict:
      #self.uc=
  def _transformZincblendeToWurtzite(self):
    az = self.uc.a
    aw = az*sqrt(2)/2. 
    cw = 2.*aw*sqrt(2./3)
    self.uc.a=aw
    self.uc.b=aw
    self.uc.c=cw 
    self.uc.gamma = 120

  def rotate_vec_parralel_to_rotax(self,vec):
    self.crystalRotMat = rotmat3Dfrom2vectors(vec,np.array([0,0,1]))
    
  def rotate_around_rotax_to_get_vecs_in_plane(self,v0,v1):
    v0 = cross(np.array([0,0,1]),v0)
    v1 = cross(np.array([0,0,1]),v1)
    #print rotmat3Dfrom2vectors(v0,v1)
    self.crystalRotMat = rotmat3Dfrom2vectors(v0,v1)*self.crystalRotMat
    
  def _get_uc_from_tuple(self,uc):
    uc_dict = dict(a=uc[0],b=uc[1],c=uc[2],alpha=uc[3],beta=uc[4],gamma=uc[5])
    uc_class = dict2class(uc_dict)
    return uc_class

  def _isallowed(self,hkl):
    if self.packing is 'fcc':
      if not isodd(sum(hkl)):
        isallowed = True
      else:
        isallowed = False
    if self.packing is 'bcc':
      if isodd(hkl[0])==isodd(hkl[1])==isodd(hkl[2]):
        isallowed = True
      else:
        isallowed = False
    if self.packing is 'diamond':
      if (isodd(hkl[0])==isodd(hkl[1])==isodd(hkl[2])) or (sum(hkl)/4.).is_integer():
        isallowed = True
      else:
        isallowed = False
    if self.packing is 'cubic':
        isallowed = True
    else:
      print "crystal structure not implemented (yet)"

  #def setCrystalRotationFromEtaRotationAxis(self):
    #"""Calculates the rotation matrix to get from the standard ctystal rotation (normal axis to uc ab plane) to the phiRotationAxis"""
    #rotnormal = cross(np.array([0,0,1]),self.etaRotationAxis)
    #rotnormal = rotnormal/numpy.linalg.norm(rotnormal)
    #rotangle = arccos(dot(np.array([0,0,1]),self.etaRotationAxis))
    #self.crystR0 = rotmat3D(rotnormal,rotangle)
  #def getRotMatFromRealHkl(self,hkl):
  #def getRotMatFromReciprocalHkl(self,hkl):

  def plotunitcell(self):
    av,bv,cv = self.ucvectors()
    av = av[0];bv=bv[0];cv=cv[0]
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    ion()
    fig = figure()
    ax = fig.gca(projection='3d')
    ax.plot([0,av[0]],[0,av[1]],[0,av[2]],'r',label='a')
    ax.plot([0,bv[0]],[0,bv[1]],[0,bv[2]],'g',label='b')
    ax.plot([0,cv[0]],[0,cv[1]],[0,cv[2]],'b',label='c')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    axis('equal')
    legend()
    draw() 

  def plotrecunitcell(self):
    av,bv,cv = self.reclattvectors()
    import matplotlib as mpl
    from mpl_toolkits.mplot3d import Axes3D
    ion()
    fig = figure()
    ax = fig.gca(projection='3d')
    ax.plot([0,av[0]],[0,av[1]],[0,av[2]],'r',label='a*')
    ax.plot([0,bv[0]],[0,bv[1]],[0,bv[2]],'g',label='b*')
    ax.plot([0,cv[0]],[0,cv[1]],[0,cv[2]],'b',label='c*')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    axis('equal')
    legend()
    draw() 

class xray(object):
  def __init__ (self,Ephot=lam2E(1),XrayDirection=[0,1,0],phiRotationAxis=[0,0,1],etaRotationAxis=[1,0,0],sampleIsVertical = False,crystal = None):
    """crystal coordinate system. The diffractometer rotation phi is assumed in crystal frame z-direction"""
    self.Ephot = Ephot
    self.Lambda = E2lam(self.Ephot)
    self.phiRotationAxis = np.array(phiRotationAxis)
    self.etaRotationAxis = np.array(etaRotationAxis)
    self.XrayDirection = np.array(XrayDirection)
    self.K = 2*pi/self.Lambda
    self.Kvec = self.getKvec()
    self._setIncidenceAngle()
    self.sampleIsVertical = sampleIsVertical
    if crystal:
      self.Xtal=crystal


  def setEphot(self,Ephot):
    self.Ephot = Ephot
    self.Lambda = E2lam(self.Ephot)
    self.K = 2*pi/self.Lambda
    self.Kvec = self.getKvec()

  def getKvec(self):
    k = 2*pi/E2lam(self.Ephot)*self.XrayDirection
    return k

  def setIncidenceAngle(self,eta=0):
    eta = eta*pi/180
    self._setIncidenceAngle(eta)

  def _setIncidenceAngle(self,eta=0):
    self.incAngrotmat = rotmat3D(self.etaRotationAxis,eta)
    self.phiRotationAxis = rotmat3D(self.etaRotationAxis,eta)*np.asmatrix(self.phiRotationAxis).ravel().transpose()
    self.phiRotationAxis = np.array(self.phiRotationAxis).ravel()
    self.eta = eta

  #def getIncidenceAngle(self,eta):
    #"""get """
    #self.etaRotationAxis = rotmat3D([0,1,0],-eta)*matrix(self.etaRotationAxis).ravel().transpose()
    #self.etaRotationAxis = np.array(self.etaRotationAxis).ravel()
    #self.eta = pi/2 - arccos(dot(self.etaRotationAxis,self.XrayDirection))

  def setPhiRotationAxis(self,v):
    self.phiRotationAxis = v/numpy.linalg.norm(v)

  def setXrayDirection(self,v):
    """Sets the X-ray direction as a vector (recommendended to be left at [0,1,0] for now)"""
    self.XrayDirection = v/numpy.linalg.norm(v)

  def getPhiRotationhkl(self,hkl,crystal=None):
    phi,phiE = self._getPhiRotationhkl(hkl,crystal=None)
    phi = phi*180/pi
    phiE = phiE*180/pi
    return phi,phiE

  def _getPhiRotationhkl(self,hkl,crystal=None):
    if not crystal:
      crystal = self.Xtal
    Q = crystal.getQhkl(hkl)
    #Q = self.crystR0*matrix(Q).ravel().transpose()
    Q = np.array(Q.ravel())
    phi,phiE = self._getPhiRotation(Q)
    return phi,phiE
#####

  def getPhiRotationhkl_new(self,hkl,crystal=None):
    phi,phiE = self._getPhiRotationhkl_new(hkl,crystal)
    phi = phi*180/pi
    phiE = phiE*180/pi
    return phi,phiE

  def _getPhiRotationhkl_new(self,hkl,crystal=None):
    Q = crystal.getQhkl(hkl)
    #Q = self.crystR0*matrix(Q).ravel().transpose()
    Q = np.array(Q.ravel())
    phi,phiE = self._getPhiRotation_new(Q)
    return phi,phiE


#####
  def getQpQn(self,Q):
    Qnorm = self.getQnorm(Q)
    Qn = dot(Q,self.phiRotationAxis)
    Qp = sqrt(Qnorm**2-Qn**2)
    return Qp,Qn
  
  def getQnQp(self,Q):
    Qnorm = self.getQnorm(Q)
    Qn = dot(Q,self.phiRotationAxis)
    Qp = sqrt(Qnorm**2-Qn**2)
    return Qn,Qp

  def getQnorm(self,Q):
    Qnorm = numpy.linalg.norm(Q)
    return Qnorm

  def getPhiRotation(self,Qcryst):
    phi,phiE =  self._getPhiRotation(Qcryst)
    return phi,phiE

  def _getPhiRotation(self,Qcryst):
    Qnorm = self.getQnorm(Qcryst)
    Qn = dot(Qcryst,np.array([0,0,1]))
    Qp = sqrt(Qnorm**2-Qn**2)
    phi = -numpy.arctan2(Qcryst[1],Qcryst[0])
    phiE = -arcsin(Qnorm**2/(2*self.K*Qp*cos(self.eta))-Qn/Qp*tan(self.eta))
    return phi,phiE

  def getPhiERotation_QnQp(self,Qn,Qp):
    phiE = self._getPhiERotation_QnQp(Qn,Qp)
    return phiE*180/pi

  def _getPhiERotation_QnQp(self,Qn,Qp):
    Qnorm = sqrt(Qp**2+Qn**2)
    phiE = -arcsin(Qnorm**2/(2*self.K*Qp*cos(self.eta))-Qn/Qp*tan(self.eta))
    return phiE

  #def _getPhiRotation_new(self,Qcryst):
    #Qnorm = self.getQnorm(Qcryst)
    #Qn = dot(Qcryst,np.array([0,0,1]))
    #Qp = sqrt(Qnorm**2-Qn**2)
    #phi = numpy.arctan2(Qcryst[1],Qcryst[0])
    #phiE = arcsin(Qnorm**2/(2*self.K*Qp*cos(self.eta))-Qn/Qp*tan(self.eta))
    #phiE = -arcsin((1/cos(self.eta)* (2* self.K* Qn* sin(self.eta)-Qp**2-Qn**2))/(2* self.K* Qp))
    #return phi,phiE

  def getDiffrationAnglesHkl(self,hkl,crystal=None):
    az,el = self._getDiffrationAnglesHkl(hkl,crystal=crystal)
    az = az*180/pi
    el = el*180/pi
    return az,el
  
  def printGeomHkl(self,hkl,crystal=None):
    phi,phiE = self.getPhiRotationhkl(hkl,crystal)
    hklstr = '[%s]' % ', '.join(map(str, hkl))
    print 'Reflection %s at phi = %f deg from preset orientation' % (hklstr,phi+phiE)
    print '  phiE = %f deg;  phi = %f deg' % (phiE,phi)

    self.getRobotAnglesHkl(hkl,crystal=None)



  def getRobotAnglesHkl(self,hkl,crystal=None):
    az,el = self._getDiffrationAnglesHkl(hkl,crystal=None)
    if self.sampleIsVertical:
      el0 = el
      az0 = az
      el,az = revertElAz(el,az)

    az = az*180/pi
    el = el*180/pi
    print 'Detector azimuth = %f deg and elevation = %f deg in Sample system' %(az0*180/pi,el0*180/pi)
    print 'Detector azimuth = %f deg and elevation = %f deg in Robot system' %(az,el)
    return az,el

  #def _getDiffrationAnglesHkl_(self,C,hkl):
    #Q = C.getQhkl(hkl)
    #phi,phiE = self._getPhiRotation(Q)
    #print phi
    #print phiE
    #Qp,Qn = C.getQpQn(Q)
    #Q = [Qn,0,Qp]
    #Q = rotmat3D(np.array([0,0,1]),phi+phiE)*matrix(Q).ravel().transpose()
    #Q = rotmat3D(self.etaRotationAxis,self.eta)*matrix(Q).ravel().transpose()
    #kout = self.Kvec+np.asarray(Q).ravel().transpose()
    #el = arcsin(kout[2]/self.K)
    #az = arctan(kout[1]/kout[0])
    #return az,el

  def _getDiffrationAnglesHkl(self,hkl,crystal=None):
    if not crystal:
      crystal = self.Xtal
    Q = crystal.getQhkl(hkl)
    phi,phiE = self._getPhiRotation(Q)
    Q = rotmat3D(np.array([0,0,1]),(phi+phiE))*np.asmatrix(Q).ravel().transpose()
    Q = rotmat3D(self.etaRotationAxis,self.eta)*np.asmatrix(Q).ravel().transpose()
    kout = self.Kvec+np.asarray(Q).ravel().transpose()
    el = arcsin(kout[2]/self.K)
    az = -1*arctan2(kout[0],kout[1])
    return az,el

def gamdel2Qfibold(gamma,delta,alpha,lam):
  gamma = np.array(iterfy(gamma))
  delta = np.array(iterfy(delta))

  shpgam = np.shape(gamma)
  shpdel = np.shape(delta)
  if not shpgam==shpdel:
    print "gamma and delta array must have same shape!"
    return
  gamma = gamma.ravel()
  delta = delta.ravel()
  Qs =  2*np.pi/lam * -np.array(rotmat3D([0,1,0],alpha)*np.mat([
    np.cos(delta)*np.cos(gamma)-1,
    -np.cos(delta)*np.sin(gamma),
    -np.sin(delta)]))
  Qip = np.sign(Qs[1,:])*np.sqrt(Qs[0,:]**2+Qs[1,:]**2);
  Qop = Qs[2,:]
  Qip = Qip.reshape(shpgam)
  Qop = Qop.reshape(shpgam)
  return Qip,Qop

def gamdel2Qfib(gamma,delta,alpha,lam):
  gamma = np.array(iterfy(gamma))
  delta = np.array(iterfy(delta))

  shpgam = np.shape(gamma)
  shpdel = np.shape(delta)
  if not shpgam==shpdel:
    print "gamma and delta array must have same shape!"
    return
  gamma = gamma.ravel()
  delta = delta.ravel()
  Qs =  2*np.pi/lam * np.array((-rotmat3D([0,1,0],-alpha))*np.mat([
    np.cos(delta)*np.cos(gamma)-1,
    -np.cos(delta)*np.sin(gamma),
    -np.sin(delta)]))
  Qip = np.sign(Qs[1,:])*np.sqrt(Qs[0,:]**2+Qs[1,:]**2);
  Qop = Qs[2,:]
  Qip = Qip.reshape(shpgam)
  Qop = Qop.reshape(shpgam)
  return Qip,Qop

def revertElAz(el,az):
  elO = np.arcsin(np.cos(el)*np.sin(az))                                                               
  azO=-np.arctan(np.sin(el)/(np.cos(el)*np.cos(az)))
  return elO,azO


