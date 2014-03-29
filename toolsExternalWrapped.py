import numpy as np
from ixppy import wrapFunc

numpyfuncs = [
    'mean',
    'std',
    'sum',
    'exp',
    'avverage',
    ]

for funName in numpyfuncs:
  exec('%s = wrapFunc(np.%s)'%(funName,funName))
