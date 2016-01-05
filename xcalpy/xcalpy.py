hasXraylib = (False,False)
try:
  import xraylib
  hasXraylib[0] = True
except:
  print 'Did not find xraylib.'
try:
  import xraylib_np
  hasXraylib[1] = True
except:
  print 'Did not find xraylib_np.'
hasPeriodictable = (False,)
try:
  import periodictable
  from periodictable import formulas 
  from periodictable.formulas import Formula as _Formula
  hasPeriodictable[0] = True
except:
  print 'Did not find periodictable.'


from scipy import constants
import numpy as np
import consts as c


class Formula(_Formula):
  def __init__(compound,density=None,natural_density=None,name=None):
    if compound == None or compound == '':
        structure = tuple()
    elif isinstance(compound,_Formula):
        structure = compound.structure
        if density is None and natural_density is None: 
            density = compound.density
        if name is None: name = compound.name
    elif _formulas.isatom(compound):
        structure = ((1,compound),)
    elif isinstance(compound,dict):
        structure = _formulas.convert_to_hill_notation(compound)
    elif _is_string_like(compound):
        try:
            formula = _formula.parse_formula(compound, table=table)
            if name: formula.name = name
            if density is not None: formula.density = density
            elif natural_density is not None: formula.natural_density = natural_density
            return formula
        except ValueError as exception:
            raise ValueError(str(exception))
            #print "parsed",compound,"as",self
    else:
        try:
            structure = _formulas._immutable(compound)
        except:
            raise ValueError("not a valid chemical formula: "+str(compound))
    _Formula.__init__(structure=structure, name=name, density=density,
                   natural_density=natural_density)
    #_Formula.__init__(*args,**kwargs)


    
def getCompoundFormula(compound,name=None,natural_density=None,density=None):
    if compound == None or compound == '':
        structure = tuple()
    elif isinstance(compound,_Formula):
        structure = compound.structure
        if density is None and natural_density is None: 
            density = compound.density
        if name is None: name = compound.name
    elif _formulas.isatom(compound):
        structure = ((1,compound),)
    elif isinstance(compound,dict):
        structure = _formulas.convert_to_hill_notation(compound)
    elif _is_string_like(compound):
        try:
            formula = _formula.parse_formula(compound, table=table)
            if name: formula.name = name
            if density is not None: formula.density = density
            elif natural_density is not None: formula.natural_density = natural_density
            return formula
        except ValueError as exception:
            raise ValueError(str(exception))
            #print "parsed",compound,"as",self
    else:
        try:
            structure = _formulas._immutable(compound)
        except:
            raise ValueError("not a valid chemical formula: "+str(compound))
    return Formula(structure=structure, name=name, density=density,
                   natural_density=natural_density)

def getElementZ(inp):
  if type(inp) is int:
    return inp
  elif type(inp) is str:
    return xraylib.SymbolToAtomicNumber(inp)

def getElementsName(name):
  return name

def attlen(material,E,density=None):
  xraylib_np.CS_Total(np.asarray(getZ(material)),np.asarray(E))

def refrIdx(material,E,density=None):
  return 

def crlAbsorption(r,sigx,sigy=None,dim=2,mu=1):#attlen('Be')):
  if sigy is None:
    sigy=sigx
  if dim==2:
    sigx_out = np.sqrt(1/(1/sigx**2+2/mu/r))
    sigy_out = np.sqrt(1/(1/sigy**2+2/mu/r))
    T = (1/sigx/sigy)*sigx_out*sigy_out
    return T,(sigx_out,sigy_out)


def crlFocalLength(r,E,material='Be',density=None):
  if density is None:
    density = xraylib.ElementDensity(xraylib.SymbolToAtomicNumber(material))
  return r/2/(1-xraylib.Refractive_Index_Re(getElementsName(material),E,density))

def crlRadius(f,E,material='Be',density=None):
  if density is None:
    density = xraylib.ElementDensity(xraylib.SymbolToAtomicNumber(material))
  #print getElementsName(material),E,density
  return f*2.*(1.-xraylib.Refractive_Index_Re(getElementsName(material),E,density))

def crlGetLensComb(r,lenses='Be',availableRs=np.asarray([50,100,200,300,500,1000,1500])*1e-6):
  availableRs.sort()
  lensset = []
  tLensIdx=0
  recr_remainder = 1./r
  while tLensIdx < len(availableRs):
    if tLensIdx==(len(availableRs)-1):
      tN = np.round(recr_remainder*availableRs[tLensIdx])
    else:
      tN = np.floor(recr_remainder*availableRs[tLensIdx])
    if tN>0:
      lensset.append((tN,availableRs[tLensIdx]))
      recr_remainder = recr_remainder - tN/availableRs[tLensIdx]
    tLensIdx+=1
  return lensset

def crlGetLensCombBinary(r,lenses='Be',availableRs=np.asarray([50,100,200,300,500,1000,1500])*1e-6):
  availableRs.sort()
  lensset = []
  tLensIdx=0
  recr_remainder = 1./r
  minrecrstep = 1./availableRs[-1]
  while tLensIdx < len(availableRs):

    if tLensIdx==(len(availableRs)-1):
      tN = np.round(recr_remainder*availableRs[tLensIdx])
      tN = int(tN>0)
      lensset.append((tN,availableRs[tLensIdx]))
    else:
      if np.abs(recr_remainder - 1./availableRs[tLensIdx]) <= minrecrstep:
	#lensset.append(True)
	lensset.append((1,availableRs[tLensIdx]))
	recr_remainder = recr_remainder - 1./availableRs[tLensIdx]
      elif recr_remainder > 1./availableRs[tLensIdx]:
	#lensset.append(True)
	lensset.append((1,availableRs[tLensIdx]))
	recr_remainder = recr_remainder - 1./availableRs[tLensIdx]
      else:
	#lensset.append(False)
	lensset.append((0,availableRs[tLensIdx]))
    tLensIdx+=1
  return lensset

def crlGetRad(lensset):
  ls = np.asarray(lensset)

  return 1./np.sum(ls[:,0]/ls[:,1])

#def getNmax(r,rphys):
    #return np.floor(recr/rphys)
#
  #getN = lambda(r)
  #recradii = np.asarray([1.*tn/tr for tr in availableRadii for tn in Nmaxeach])
#
  #recradii.sort()
#
  #availableRadii


  

