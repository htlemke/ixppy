import numpy as np
from ixppy import wrapFunc

numpyfuncs = [
    'mean',
    'std',
    'sum',
    'exp',
    'nansum',
    'average',
    'apply_over_axes'
    ]

for funName in numpyfuncs:
  exec('%s = wrapFunc(np.%s)'%(funName,funName))

def Nansum(data):
  o = np.squeeze(np.apply_over_axes(np.nansum,data,range(1,np.rank(data))))
  return np.atleast_2d(o).T

nansum = wrapFunc(Nansum)
