import numpy as np
def rebin(a, shape):
    sh = shape[0],a.shape[0]//shape[0],shape[1],a.shape[1]//shape[1]
    return a.reshape(sh).mean(-1).mean(1)

def edgeMask(imgOrShape,nPixels):
  """ mask the edge of the images """
  if (isinstance(imgOrShape,list) or isinstance(imgOrShape,tuple)):
    mask = np.ones(imgOrShape,dtype=np.bool)
  else:
    mask = np.ones_like(imgOrShape,dtype=np.bool)
  if ndim == 1:
    mask[:nPixels]  = False
    mask[-nPixels:] = False
  elif ndim == 2:
    mask[:nPixels,:]  = False
    mask[:,:nPixels]  = False
    mask[:,-nPixels:] = False
    mask[-nPixels:,:] = False
  elif ndim == 3:
    for m in mask:
      m[:nPixels,:]  = False
      m[:,:nPixels]  = False
      m[:,-nPixels:] = False
      m[-nPixels:,:] = False
  else:
    print "Cannot make mask"
  return mask
