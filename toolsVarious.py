import types
import os
import numpy as np
from copy import deepcopy
import re

def fileExists(fname):
	return os.path.exists(fname)

def removeFileIfExists(fname):
	if (fileExists(fname)):
		os.remove(fname)

def h5group2class(d):
  import h5py
  if (not isinstance(d,h5py.Group) ):
    return None
  else:
    class Dummy(object):
      pass  
    c = Dummy()
    for elem in d.keys():
        c.__dict__[elem] = d[elem].value
    return c
   

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


def isodd(num):
  return num & 1 and True or False


def iterfy(iterable):
    if isinstance(iterable, basestring):
        iterable = [iterable]
    try:
        iter(iterable)
    except TypeError:
        iterable = [iterable]
    return iterable

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

def sind(x): return np.sin(np.rad2deg(x))
def cosd(x): return np.cos(np.rad2deg(x))
def tand(x): return np.tan(np.rad2deg(x))
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


def dict_merge(a, b):
    '''recursively merges dict's. not just simple a['key'] = b['key'], if
    both a and bhave a key who's value is a dict then dict_merge is called
    on both values and the result stored in the returned dictionary.'''
    if not isinstance(b, dict):
        return b
    result = deepcopy(a)
    for k, v in b.iteritems():
        if k in result and isinstance(result[k], dict):
                result[k] = dict_merge(result[k], v)
        else:
            result[k] = deepcopy(v)
    return result

def varName(varStr): 
  return re.sub('\W|^(?=\d)','_', varStr)
