import numpy as np
from ixppy import wrapFunc
import copy

numpyfuncs = [
    'mean',
    'std',
    'sum',
    'exp',
    'nansum',
    'average',
    'apply_over_axes',
    'polyval',
    'isnan',
    'any'
    ]

for funName in numpyfuncs:
  exec('%s = wrapFunc(np.%s)'%(funName,funName))

def Nansum(data):
  o = np.squeeze(np.apply_over_axes(np.nansum,data,range(1,np.rank(data))))
  return np.atleast_2d(o).T


nansum = wrapFunc(Nansum,transposeStack=False)


def MaskNan(data,mask):
  if len(data)>0:
    dat = data.astype(np.float)
    dat[:,mask] = np.nan
    return dat
  else:
    return data
maskNan = wrapFunc(MaskNan)

