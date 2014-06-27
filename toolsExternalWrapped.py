import numpy as np
from ixppy import wrapFunc

numpyfuncs = [
    'mean',
    'std',
    'sum',
    'exp',
    'average',
    'apply_over_axes'
    ]

for funName in numpyfuncs:
  exec('%s = wrapFunc(np.%s)'%(funName,funName))
