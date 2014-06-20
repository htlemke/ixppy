import types
import os
import numpy as np
from copy import deepcopy
import re
from itertools import takewhile
from toolsVecAndMat import smooth
#from ixppy import wrapFunc

def allnamesequal(name):
  return all(n==name[0] for n in name[1:])

def commonPathPrefix(paths, sep='/'):
  bydirectorylevels = zip(*[p.split(sep) for p in paths])
  return sep.join(x[0] for x in takewhile(allnamesequal, bydirectorylevels))

class dropObject(object):
  def __init__(self,name='noname',parent=None):
    self._name = name
    self._parent = parent
    self._ixpsaved = []

  def _add(self,name,data,ixpsaved='auto'):
    self.__dict__[name]=data
    self._ixpsaved.append((name,ixpsaved))
  def __repr__(self):
    return "dropObject with fields: "+str(self.__dict__.keys())
  def __getitem__(self,x):
    return self.__dict__[x]
  def __setitem__(self,name,var,setParent=True):
    self._add(name,var)
    if setParent:
      try:
        self[name]._parent = self
      except:
	pass
  def _get_keys(self):
    return [tk for tk in self.__dict__.keys() if not tk[0]=='_']

def itemgetToIndices(x,size,boolean=False):
  if size==0:
    return []
  if (type(x) is int) or (type(x) is np.int64):
    xo = np.array([x])
  elif (type(x) is tuple) and isinstance(x[0],np.ndarray):
    xo = x[0]
  elif isinstance(x,slice):
    xo = np.arange(*x.indices(size))
  else:
    raise IndexError
    
  if (max(xo)>=size):
    raise IndexError
  if boolean:
    xo = np.zeros(size,dtype=bool)[xo] = True
  return xo

def getFromObj(obj,name):
  temp = name
  where = obj
  while (temp.find(".")>0):
    parent = temp[0:temp.find(".")]
    where =  where.__dict__[parent]
    temp = temp[temp.find(".")+1:]
  return where.__dict__[temp]

def existsInObj(obj,name):
  temp = name
  where = obj
  while (temp.find(".")>0):
    parent = temp[0:temp.find(".")]
    if (parent in where.__dict__):
      where =  where.__dict__[parent]
      temp = temp[temp.find(".")+1:]
    else:
      return False
  return (temp in where.__dict__)


def addToObj(obj,name,value,overWrite=True,ixpsaved=False, setParent=True):
  """ Functions to add things to an object create intermediate dropObject if 
  necessary
  usage:
  from a class you could call addToObj(self,"p1.p2.p3.a1",np.arange(10))
  p1 p2 and p3 will be dropObject and a1 an attribute of p3 with value np.arange(10)
  given an instance (for example data) of a class another use is addToObj(data,...)
  """
  temp = name
  where = obj
  while (temp.find(".")>0):
    parent = temp[0:temp.find(".")]
    if parent not in where.__dict__:
      where.__dict__[parent] = dropObject(name=parent,parent=where)
    if ixpsaved:
      if not hasattr(where,'_ixpsaved'):
	where._ixpsaved = []
      if not (parent in [ix[0] for ix in where._ixpsaved]):
        where._ixpsaved.append((parent,'auto'))
    where =  where.__dict__[parent]
    temp = temp[temp.find(".")+1:]
  if ( (temp not in where.__dict__) or overWrite ):
    if hasattr(where,'_ixpsaved'):
      where[temp] = value
    else:
      where.__dict__[temp] = value
    if setParent:
      try:
        where.__dict__[temp]._parent = where
      except:
	pass
  return obj

def fileExists(fname):
	return os.path.exists(fname)

def removeFileIfExists(fname):
	if (fileExists(fname)):
		os.remove(fname)

def h5GroupToObj(d):
  import h5py
  if (not isinstance(d,h5py.Group) ):
    return None
  else:
    c = dropObject()
    for elem in d.keys():
        c._add(elem,d[elem][...])
    return c

def isIterable(something):
  return hasattr(something,'__iter__')

def iterate(data,function,*args,**keywords):
  if (not isIterable(data)):
    data=iterfy(data)
  nargs = len(args)
  nkey  = len(keywords)
  if   ( (nargs!=0) and (nkey!=0) ):
    return [function(x,args,keywords) for x in data]
  elif ( (nargs!=0) and (nkey==0) ):
    return [function(x,*args) for x in data]
  elif ( (nargs==0) and (nkey!=0) ):
    return [function(x,**keywords) for x in data]
  else:
    return [function(x) for x in data]

def iterdepth(iterable):
  """only for lists/arrays, only along first element"""
  if isiter(iterable):
    N = 0
    iter = True
    while iter:
      N+=1
      iter = eval('isiter(iterable'+(N)*'[0]'+')')
  else: 
    N=0
  return N

def strucArrayToObj(data):
  """
  Transform a structured array as class
  x = np.zeros(3, dtype=[('x','f4'),('y',np.float32),('value','f4',(2,2))])
  A=strucArrayToObj(x)
  print A.value
  """
  c = dropObject()
  if data[0].dtype.names is not None:
    for fieldname in data[0].dtype.names:
      c._add(fieldname,data[fieldname])
  else:
    print "No clue on how to make an object out handle ",data
  return c

def dictToObj(d):
    """Return a class that has same attributes/values and 
       dictionaries key/value
    """
    #see if it is indeed a dictionary
    if type(d) != types.DictType:
        return None
    c = dropObject()
    for elem in d.keys():
        c._add(elem,d[elem])
    return c

def isodd(num):
  return num & 1 and True or False

def get_copies(obj):                                                            
    return [objname for objname,oid in globals().items() if id(oid)==id(obj)]  

def iterfy(iterable):
    if isinstance(iterable, basestring):
        iterable = [iterable]
    try:
        iter(iterable)
    except TypeError:
        iterable = [iterable]
    return iterable

def iterDiff(it1,it2):
  if (it1 is None) or (it2 is None):
    diff = set([])
  else:
    diff = (set(it1)-set(it2))
  return diff

def num2sci(num,unit='',precision=3):
  exponents = np.arange(-24,24+1,3)
  prefixes =['y','z','a','f','p','n','u','m','','k','M','G','T','P','E','Z','Y']
  exponNum = np.floor(np.log10(num))
  exponSci = np.floor(exponNum/3)*3
  if not exponSci in exponents:
    exponSci = exponents[np.argmin(np.abs(exponents-exponSci))]
  sci = num / 10**exponSci
  prefix = prefixes[list(exponents).index(exponSci)]
  sci = round(sci,precision)
  if unit is '':
    return sci,exponSci
  else:
    return str(sci)+prefix+unit

asin = np.arcsin
acos = np.arccos
atan = np.arctan
atan2 = np.arctan2

def sind(x): return np.sin(np.deg2rad(x))
def cosd(x): return np.cos(np.deg2rad(x))
def tand(x): return np.tan(np.deg2rad(x))
def arcsind(x): return np.rad2deg(np.arcsin(x))
def arccosd(x): return np.rad2deg(np.arccos(x))
def arctand(x): return np.rad2deg(np.arctan(x))
def arctan2d(x): return np.rad2deg(np.arctan2(x))
asind = arcsind
acosd = arccosd
atand = arctand
atan2d = arctan2d


def round_err(value,err,asstring=False):
	""" returns value and err rounded to the number of
			significant digits of err """
	if (not (np.isfinite(err))):
		err = np.abs(value/1e3)
	if (err != 0):
		ndigerr = -int(np.floor(np.log10(err)))
		if (ndigerr<1): ndigerr=2
		v =np.around(value,ndigerr)
		e =np.around(err,ndigerr)
	else:
		v=value
		e=err
	if (asstring):
		return "%s" % v,"%s" % e
	else:
		return v,e

def saveTXT(fname,x,Ys,form="%+.6g",sep=" ",header=None,headerv=None):
	""" Write data to file 'fname' in text format.
			Inputs:
				x = x vector
				Ys = vector(or array or list of vectors) for the Ys
				form = format to use
				sep = separator
				header = text header (must be a string)
				headerv = vector to be used as header, it is convienient when
					the output must be of the form
						Ncol 252 253 254
						x1   y11 y12 y13
						.......
					In this case headerv should be [252,253,254]
	"""
	if (type(x) != np.ndarray): x=np.array(x)
	if (type(Ys) != np.ndarray): Ys=np.array(Ys)
	if (len(Ys.shape)==1):
		Ys=Ys.reshape(Ys.shape[0],1)
	nx = len(x)
	if (Ys.shape[0] == nx):
		ny=Ys.shape[1]
	elif (Ys.shape[1] == nx):
		ny=Ys.shape[0]
		Ys=np.transpose(Ys)
	else:
		raise Exception("dimension of x (%d) does not match any of the dimensions of Ys (%d,%d)" % (nx,Ys.shape[0],Ys.shape[1]))
	try:
		import codecs
		f=codecs.open(fname,encoding='utf-8',mode="w")
	except ImportError:
		f=open(fname,"w")
	if (header is not None):
		f.write(header.strip()+"\n")
	if (headerv is not None):
		f.write("%d" % (ny+1))
		for i in range(ny):
			f.write(sep)
			f.write(form % headerv[i])
		f.write("\n")
	for i in range(nx):
		f.write(form % x[i])
		f.write(sep)
		for j in range(ny-1):
			f.write(form % Ys[i,j])
			f.write(sep)
		f.write(form % Ys[i,-1])
		f.write("\n")
	f.close()


def dictMerge(a, b):
    '''recursively merges dict's. not just simple a['key'] = b['key'], if
    both a and bhave a key who's value is a dict then dict_merge is called
    on both values and the result stored in the returned dictionary.'''
    if not isinstance(b, dict):
        return b
    result = deepcopy(a)
    for k, v in b.iteritems():
        if k in result and isinstance(result[k], dict):
                result[k] = dictMerge(result[k], v)
        else:
            result[k] = deepcopy(v)
    return result

def varName(varStr): 
  return re.sub('\W|^(?=\d)','_', varStr)

def find_peaks_zc(x,smoothWindow=10,min_dist=None,max_dist=None):
  xd = smooth(np.diff(x),smoothWindow)
  zc = np.where(np.diff(np.sign(xd)))[0]
  Dzc = np.diff(zc)
  if min_dist is not None:
    sel = Dzc>min_dist
    Dzc = Dzc[sel]
    selind = (~sel).nonzero()[0]
    zc = np.delete(zc,np.hstack([selind,selind+1]))
  if max_dist is not None:
    sel = Dzc<max_dist
    Dzc = Dzc[sel]
    selind = (~sel).nonzero()[0]
    zc = np.delete(zc,np.hstack([selind,selind+1]))
  return zc,Dzc 





#polyval = wrapFunc(np.polyval)
def dict2class(d):
    """Return a class that has same attributes/values and 
       dictionaries key/value
    """
    
    #see if it is indeed a dictionary
    if type(d) != types.DictType:
        return None
    
    #define a dummy class
    class Dummy:
        pass
        
    c = Dummy
    for elem in d.keys():
        c.__dict__[elem] = d[elem]
    return c

def corrNonlinGetPar(data,correct,order=2,data_0=0,correct_0=0):
  p =  np.polyfit(data,correct,order)
  p[-1] = p[-1]-correct_0
  return p

def corrNonlin(data,polypar,data_0=0,correct_0=0):
  m = 1/np.polyval(np.polyder(polypar),data_0)
  return m*(np.polyval(polypar,data)-correct_0) + data_0


