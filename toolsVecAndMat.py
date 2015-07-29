import numpy as np
import numpy.linalg as linalg
import time
import types
import numpy.ma as ma
from toolsLog import logbook
from scipy import percentile
import toolsDistrAndHist,tools
import copy

def smartIdx(idx,forceContigous=False):
  """ Try to interpret an array of bool as slice;
  this allows selecting a subarray alot more efficient 
  since array[slice] it returns a view and not a copy """
  if (isinstance(idx,int)):
    ret = slice(idx,idx+1)
  else:
    idx = np.asarray(idx)
    if idx.dtype == np.bool:
      i = np.where(idx)[0]
    else:
      i = idx
    # in case there is only one
    if (len(i) == 1):
      ret = slice(i[0],i[0]+1)
      return ret
    if forceContigous:
      ret = slice(i[0],i[-1])
    else:
      d = i[1:]-i[0:-1]
      dmean = int(d.mean())
      if np.all(d==dmean):
        ret = slice(i[0],i[-1]+1,dmean)
      else:
        ret = idx
  return ret

def filtvec(v,lims,getbool=False):
  if not getbool:
    return np.logical_and(v>min(lims),v<max(lims))
  else:
    return  (v>min(lims),v<max(lims))

def subset(a,lims):
  """cut a region from 2D array a given by ginput lims"""
  lims = np.round(lims)
  if len(np.shape(lims))==1 and len(lims)==4:
    lims = np.reshape(lims,[2,2]).T

  return a[min(lims[:,1]):max(lims[:,1]) , min(lims[:,0]):max(lims[:,0])]

def filtWithPercentile(v,low,high):
	""" return indeces that satisty the property that
	values for this indeces are within `low` and `high` in the probability 
	distribution """
	lims = percentile(v, (low,high) )
	return filtvec(v,lims)

def filterWithMad(v,fac=1,mad=None):
	""" returns boolean indeces that are true when they satify the
	relation of 'v' being within its median +/- 'fac' times the RMS calculated
	with the MAD """
	if (mad is None): mad = toolsDistrAndHist.mad(v)
	return (np.abs(v-np.median(v))<fac*mad)

class AverageProfile(object):
  def __init__(self,xs,xe,dx):
    self.xs = xs
    self.xe = xe
    self.dx = dx
    self.xout = np.arange(xs+dx/2.,xe+dx/2.,dx)
    self.n  = len(self.xout)
    self.y = []
    self.ncurves = 0
    for i in range(self.n):
      self.y.append([])

  def add(self,xdata,ydata,sigma=None):
    """ sigma (if given are use to weight the average) """
    self.ncurves += 1
    for i in range(len(xdata)):
      bin = int( (xdata[i]-self.xs)/self.dx )
      if ((bin>=0) & (bin<len(self.y)-1)):
        if (sigma is None):
          self.y[bin].append( (ydata[i],1) )
        else:
          self.y[bin].append( (ydata[i],sigma[i]) )

  def get(self):
    yout = np.empty_like(self.xout)
    sout = np.empty_like(self.xout)
    nout = np.empty_like(self.xout,dtype=np.int)
    for i in range(self.n):
      # print i,self.n,self.y[i]
      if (len(self.y[i])==0):
        yout[i] = np.nan 
        sout[i] = np.nan 
        nout[i] = 0 
      else:
        temp = np.asarray(self.y[i])
				# temp[:,0] are Ys, temp[:,1] are sigmas
        yout[i] = np.sum(temp[:,0]/(temp[:,1]**2))/np.sum(1./(temp[:,1]**2)) 
        nout[i] = len(self.y[i]) 
        sout[i] = np.sqrt(1./np.sum(1./(temp[:,1]**2)))
        #sout[i] /= np.sqrt(nout[i])
    return (self.xout,yout,sout,nout)

class smartAverage(object):
  """ Average several curves (xi,yi) where the xi might be different """
  def __init__(self,xbins_edges):
    self.xedges = xbins_edges
    self.x    = (self.xedges[1:] + self.xedges[:-1])/2.
    self._w   = np.zeros_like(self.x)
    self._n   = np.zeros_like(self.x,dtype=np.uint)
    self.y    = np.zeros_like(self.x)
    self.err  = np.zeros_like(self.x)*np.nan
    self.ncurves = 0

  def add(self,x,y,sigma=None):
    """ sigma (if given are use to weight the average, as 1/sig) """
    self.ncurves += 1
    y = y.astype(np.float)
    idx = np.digitize(x,self.xedges)-1
    if sigma is None: sigma = np.ones_like(x)
    for i,bin in enumerate(idx):
      if (bin>=len(self.y)): continue
      w = 1./sigma[i]**2
      self.y[bin]  = (self.y[bin]*self._w[bin] + w*y[i])/ (self._w[bin]+w)
      self._w[bin] += w
      self._n[bin] += 1
    idx = self._n>=1
    self.err[idx] = np.sqrt( 1./self._w[idx] )

  def getValidIdx(self):
    """ return index of booleans, True for bins with at least one value """
    return self._n>=1

  def getFiniteIdx(self):
    """ return index of booleans, True for bins with at least one value and finite err and y"""
    return self.getValidIdx() & np.isfinite(self.y) & np.isfinite(self.err)


def rotmat3D(v,ang):
  """3D rotation matrix around axis v about angle ang"""
  ux = v[0]
  uy = v[1]
  uz = v[2]
  c = np.cos(ang)
  s = np.sin(ang)
  rotmat = np.matrix(
  [[ux**2+(1-ux**2)*c , ux*uy*(1-c)-uz*s , ux*uz*(1-c)+uy*s],
  [ux*uy*(1-c)+uz*s , uy**2+(1-uy**2)*c , uy*uz*(1-c)-ux*s],
  [ux*uz*(1-c)-uy*s , uy*uz*(1-c)+ux*s , uz**2+(1-uz**2)*c]]);
  rotmat = np.matrix(rotmat)
  return rotmat

#def rotmat3Dfrom2vectors(v0,v1):
  #"""calculate 3D rotation matrix that rotates from v0 to v1"""
  #v0 = v0/norm(v0)
  #v1 = v1/norm(v1)
  #ax = cross(v0,v1);
  #ang = arcsin(norm(ax))
  #ax = ax/norm(ax)
  #rotmat = rotmat3D(ax,ang)
  #return rotmat

def rotmat3Dfrom2vectors(v0,v1):
  """calculate 3D rotation matrix that rotates from v0 to v1"""
  v0 = v0/linalg.norm(v0)
  v1 = v1/linalg.norm(v1)
  ax = linalg.cross(v0,v1)
  if not linalg.norm(ax)==0.:
    ax = ax/linalg.norm(ax)
    ve = linalg.cross(ax,v0)
    cx = linalg.dot(v1,v0)
    cy = linalg.dot(v1,ve)
    ang = np.arctan2(cy,cx)
    rotmat = rotmat3D(ax,ang)
  else:
    rotmat = np.eye(3)
  return rotmat

def pol2cart(theta, radius, units='rad'):
    """Convert from polar to cartesian coordinates 
     
    **usage**: 
        x,y = pol2cart(theta, radius, units='deg') 
    """
    if units in ['deg', 'degs']:
        theta = theta*np.pi/180.0
    xx = radius*np.cos(theta)
    yy = radius*np.sin(theta)
    return xx,yy
#---------------------------------------------------------------------- 
def cart2pol(x,y, units='rad'):
    """Convert from cartesian to polar coordinates 
     
    **usage**: 
        theta, radius = pol2cart(x, y, units='deg') 
         
    units refers to the units (rad or deg) for theta that should be returned"""
    radius= np.hypot(x,y)
    theta= np.arctan2(y,x)
    if units in ['deg', 'degs']:
        theta=theta*180/np.pi
    return theta, radius


def oversample(v,fac=None,interval=None):
  assert (fac is not None) or (interval is not None), "You need to define either oversampling factor or an interval"
  if fac is not None:
    vo = np.linspace(min(v),max(v),v.shape[0]*fac)
  if interval is not None:
    vo = np.arange(min(v),max(v),interval)
  return vo


class ArrayWithMasks(np.ndarray):
	""" to try it:
	a=np.arange(100).reshape(20,5)
	b=ArrayWithMasks(a)
	b.addMask("m1",b>19)
	b.addMask("m2",b<60)
	print b.m1; # masked using only mask `m1`
	print b.m2;
	print b.mdata; # masked with all masks ..
	"""
	def __new__(cls, input_array):
    # Input array is an already formed ndarray instance
    # We first cast to be our class type
		obj = np.asarray(input_array).view(cls)
    # add the new attribute to the created instance
    # Finally, we must return the newly created object:
		return obj 

	def __init__(self,a):
		self.mdata  = a
		self.mask   = np.ones(a.shape,dtype=np.bool)
		self._masks = {}
		self._infos = {}
		self._fracs = {}


	def addMask(self,name,mask,info=None):
		self._masks[name]=mask
		self._infos[name]=info
		self._fracs[name]=np.sum(mask)/float(self.size)
		self.__dict__["mdata_%s"%name] = self[mask]
		self.mask = self.mask & mask; # update total mask
		self.mdata = self[self.mask]

	def getArray(self,masks="all"):
		if (masks is None):
			return self
		elif (masks == "all"):
			return self.mdata
		elif (masks in self._masks):
			return self.__dict__[masks]
		else:
			mask = np.ones(self.shape,dtype=np.bool)
			for m in masks:
				mask = mask & self._masks[m]
			return self[mask]

	def getMaskNames(self):
		return self._masks.keys()

	def getMask(self,mask="all"):
		if (mask =="all"):
			return self.mask
		elif (mask in self._masks):
			return self._masks[mask]
		else:
			logbook("mask %s not present, returning None" % mask)
			return None

	def getInfo(self):
		s = ""
		for m in self.getMaskNames():
			if ( self._infos[m] is not None ):
				s+="# mask %s, true for %s, %s\n" % (m,self._fracs[m],self._infos[m])
			else:
				s+="# mask %s, true for %s\n" % (m,self._fracs[m])
		return s[0:-1]

def rollingFunction(v,fun,ptrad=10,truncate=False):
  v = np.array(v)
  lv = len(v)
  if not truncate:
    return np.array([fun(v[max(n-ptrad,0):min(n+ptrad+1,lv+1)]) for n in range(lv)])
  else:
    return np.array([fun(v[n-ptrad:n+ptrad+1]) for n in range(ptrad,lv-ptrad+1)])

def rollingAverage(a, n=3):
  ret = np.cumsum(a, dtype=float)
  return (ret[n-1:] - ret[:1 - n]) / n

def smooth(x,window_len=11,window='hanning'):
        if x.ndim != 1:
                raise ValueError, "smooth only accepts 1 dimension arrays."
        if x.size < window_len:
                raise ValueError, "Input vector needs to be bigger than window size."
        if window_len<3:
                return x
        if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
                raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"
        s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
        if window == 'flat': #moving average
                w=np.ones(window_len,'d')
        else:  
                w=eval('np.'+window+'(window_len)')
        y=np.convolve(w/w.sum(),s,mode='same')
        return y[window_len:-window_len+1]

def ndmesh(*args):
  args = map(np.asarray,args)
  return np.broadcast_arrays(*[x[(slice(None),)+(None,)*i] for i, x in enumerate(args)])


def flatten(x):
  if np.iterable(x):
    return [a for i in x for a in flatten(i)]
  else:
    return [x]


def rebin(a, *args):
    '''rebin ndarray data into a smaller ndarray of the same rank whose dimensions
    are factors of the original dimensions. eg. An array with 6 columns and 4 rows
    can be reduced to have 6,3,2 or 1 columns and 4,2 or 1 rows.
    example usages:
    >>> a=rand(6,4); b=rebin(a,3,2)
    >>> a=rand(6); b=rebin(a,2)
    '''
    shape = a.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(args)
    evList = ['a.reshape('] + \
             ['args[%d],factor[%d],'%(i,i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in range(lenShape)] + \
             ['/factor[%d]'%i for i in range(lenShape)]
    print ''.join(evList)
    return eval(''.join(evList))

def polyVal(comps,i0):                                                          
  i0 = np.asarray(tools.iterfy(i0))                                             
  pol = np.vander(i0,len(comps))                                        
  return np.asarray(np.matrix(pol)*np.matrix(comps.reshape((len(comps),-1)))).reshape((len(i0),)+np.shape(comps)[1:])                            
                                                                                
def polyDer(comps,m=1):                                                         
  compsf = comps.reshape((len(comps),-1))
  n = len(compsf) - 1                                                            
  y = compsf.reshape((n+1,-1))[:-1] * np.expand_dims(np.arange(n, 0, -1),1)                      
  if m == 0:                                                                    
    val = comps
    return val
  else:                                                                         
    val = polyDer(y, m - 1) 
    return val.reshape((n,)+np.shape(comps)[1:])

def polyFit(i0,Imat,order=3, removeOrders=[]):
  Imatf = Imat.reshape((len(Imat),-1))
  pol = np.vander(i0,order+1)                                                   
  removeOrders = tools.iterfy(removeOrders)                                     
  removeOrders = np.sort(removeOrders)[-1::-1]                                  
  for remo in removeOrders:                                                     
    #print remo                                                                  
    pol = np.delete(pol,-(remo+1),axis=1)                                       
  lhs = copy.copy(pol)                                                          
  scale = np.sqrt((lhs*lhs).sum(axis=0))                                        
  lhs /= scale                                                                  
  comps,resid,rnk,singv = linalg.lstsq(lhs,Imatf)                                
  comps = (comps.T/scale).T                                                     
                                                                                
  for remo in removeOrders:                                                     
    comps = np.insert(comps,order-remo,0,axis=0)                                
  return comps.reshape((order+1,)+np.shape(Imat)[1:])


def unmaskData(data,mask,fillValue=np.nan):
  """
  puts data back into its original length gven a mask of that original length.
  mask must have as many False values as the length of data.
  Works on stacks of data. Masked entries (marked by True in mask) will be 
  filled with fillValue (default is np.nan).
  """
  re = np.ones((len(data),len(mask)))*fillValue
  msk = mask.nonzero()
  for n,dat in enumerate(data):
    re[n].put(msk,dat)
  #return re.reshape([len(data),32,185,388])
  return re

def applyFuncOnClosestN(func,N,x0,x1,D0):
  acni = _ApplyClosestN(func,N,x0,x1,D0)
  return acni

def rollaxis2end(a):
  n = a.ndim
  axes = list(range(0, n))
  axes.remove(0)
  axes.insert(n, 0)
  return a.transpose(axes)

def rollaxis2beg(a):
  n = np.ndim(a)
  axes = list(range(0, n))
  axes.remove(n-1)
  axes.insert(0, n-1)
  return a.transpose(axes)


class _ApplyClosestN(object):
  def __init__(self,func,N,x0,x1,D0):
    self.x0 = np.sort(x0)
    self.x0_sortind = np.argsort(x0)
    self.D0 = D0
    self.x1 = np.sort(x1)
    self.x1_sortind = np.argsort(x1)
    self.func = func
    self.N = N
    self.presInd = range(N)
    self.presLst = list(self.x0[:N])
    self.res = []
  
  def step(self,x):
    while (self.presLst[-1]-x) < (x-self.presLst[0]):
      if self.presLst[-1]>=len(self.x0):
	break
      if (self.x0[self.presInd[-1]+1] - x) > (x-self.presLst[0]):
        self.res.append(self.func(self.D0[self.x0_sortind[np.ix(self.presInd)]]))
        return
      self.presInd.pop(0)
      self.presInd.append(self.presInd[-1]+1)
      self.presLst.pop(0)
      self.presLst.append(self.x0[self.presInd[-1]])

  def run(self):
    for x in self.x1:
      self.step(x)

def arrayFromList(a):
  shps = [np.shape(ta) for ta in a]
  shplens = np.asarray([len(tshp) for tshp in shps])
  shp = shps[shplens.argmax()]
  o = np.nan*np.ones((len(a),)+shp)
  for n,ta in enumerate(a):
    try:
      o[n] = ta
    except:
      pass
  return o  















