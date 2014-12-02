import numpy as np
from tools import getEdges

def anaPumpProbePar(I,probeMonitor,timeDelay=None,pumpMonitor=None,xrayon=None,laseron=None,delayBin=None,normPump=False):
  """General tool to analyze pump probe data in ixppy. It takes 
  care of normalization,off shots and pump intensity monitors. The result 
  is provided as pumped to unpumped ratio.
  """

  xrayonFilter   = I.ones()
  laseronFilter  = I.ones()
  laseroffFilter = I.ones()
  if xrayon is not None: xrayonFilter *= xrayon.filter(True).ones()
  if laseron is not None: 
    laseronFilter *= laseron.filter(True).ones()
    laseroffFilter *= laseron.filter(False).ones()
  filtXonLon     = xrayonFilter*laseronFilter
  filtXonLoff    = xrayonFilter*laseroffFilter
  ratioOn        = I/probeMonitor*filtXonLon
  if laseron is not None:
    ratioOff       = I/probeMonitor*filtXonLoff
    diffRatio      = ratioOn/ratioOff.median()
  else:
    diffRatio = ratioOn
  if normPump:
    pumpMonitorFiltered =  filtXonLon*pumpMonitor *I.ones()*probeMonitor.ones()
    diffRatio = (diffRatio-1)/pumpMonitorFiltered * np.mean(pumpMonitorFiltered.median())+1
  if timeDelay is not None:
    if delayBin is not None:
      if type(delayBin) is np.ndarray:
	edges = delayBin
      elif type(delayBin) is dict:
	edges = getEdges(timeDelay.R,delayBin)
      td = timeDelay.digitize(edges)
    else:
      td = timeDelay
    res = td.ones()*diffRatio
  else:
    res = diffRatio
  return res
