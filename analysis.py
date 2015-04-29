import numpy as np
from tools import getEdges,histVec
from toolsExternalWrapped import nansum
import ixppy
from nonlinearParameterCorrection import CorrNonLin

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


class PumpProbeNDarray(object):
  def __init__(self,
      data, monitor=None, timeDelay=None,
      pumpOn=None,probeOn=None,
      timeBinSize=20e-15,
      attachDropObject=None,
      resFina='tempFile_pumpProbeND.ixp.h5',
      dsNamePrefix='pumpProbe'):

    self.data = data
    self._monitor = monitor
    self._timeDelay = timeDelay
    self._pumpOn = pumpOn
    self._probeOn = probeOn
    self.timeBinSize = timeBinSize
    self._res = attachDropObject
    self.resultFile = resFina
    self._dsNamePrefix = dsNamePrefix

  def _getAttachDropObj(self):
    if self._res is None:
      self._res = ixppy.dataset(self.resultFile)
      self._res['det'] = ixppy.tools.dropObject()
    return self._res
  res = property(_getAttachDropObj)

  def _getMonitor(self): 
    if self._monitor is None:
      try:
        self._monitor = self.res[self._dsNamePrefix+'_nansum_i0']
      except:
        self.res[self._dsNamePrefix+'_nansum_raw'] = nansum(self.data)
        print "Extracting I0 from radial profile..."
        self.res[self._dsNamePrefix+'_nansum_raw'].evaluate()
        self.res[self._dsNamePrefix+'_nansum_i0'] = \
	    self.res[self._dsNamePrefix+'_nansum_raw'].get_memdata()[0]
        self._monitor = self.res[self._dsNamePrefix+'_nansum_i0']
    return self._monitor
  def _setMonitor(self,monitor):
    self._monitor = monitor
  monitor = property(_getMonitor,_setMonitor)
  
  def _getPumpOn(self):
    return self._pumpOn.filter(True).ones()
  pumpOn = property(_getPumpOn)
  def _getPumpOff(self):
    return self._pumpOn.filter(False).ones()
  pumpOff = property(_getPumpOff)
  def _getProbeOn(self):
    return self._probeOn.filter(True).ones()
  probeOn = property(_getProbeOn)
  def _getProbeOff(self):
    return self._probeOn.filter(False).ones()
  probeOff = property(_getProbeOff)

  def _getTimeDelay(self):  
    pumpedAndProbed = self.probeOn*self.pumpOn
    #pumpedAndProbed = self.probeOn*self.pumpOn*self.monitor.ones()
    if self._timeDelay is None:
      try:
        self._timeDelay = self.data.ones()*self.data.scan[0]
        timeDelayEdges = histVec(self.data.scan[0])*pumpedAndProbed
      except:
        "did not find scan variable, will use integers instead"
        self._timeDelay = self.data.ones()*range(len(self.data))
        timeDelayEdges = histVec(range(len(self.data))) 
      td = self._timeDelay.digitize(timeDelayEdges)*pumpedAndProbed
    else:
      if self.timeBinSize is not None:
        td = self._timeDelay.digitize(self.timeBinSize)*pumpedAndProbed
      else:
	td = self._timeDelay*pumpedAndProbed
    self.res[self._dsNamePrefix+'_timeDelay_binned'] = td
    return td
  timeDelay = property(_getTimeDelay)

  def _getGooddata(self):
    gooddata = self.timeDelay.ones()*self.monitor.ones()
    if xrayon is not None:
      gooddata *= xrayon.filter(True).ones()
    return gooddata
  isGoodData = property(_getGooddata)

  def _getDataOffNorm(self):
    try:
      return self.res[self._dsNamePrefix+'_UnPumpedNorm']
    except:
      self.res[self._dsNamePrefix+'_UnPumpedNorm'] = self.data*self.pumpOff*self.probeOn/self.monitor
      return self.res[self._dsNamePrefix+'_UnPumpedNorm']
  dataOffNorm = property(_getDataOffNorm)
    
  def _getDataOff(self):
    try:
      return self.res[self._dsNamePrefix+'_UnPumped']
    except:
      self.res[self._dsNamePrefix+'_UnPumped'] = self.data*self.pumpOff*self.probeOn
      return self.res[self._dsNamePrefix+'_UnPumped']
  dataOff = property(_getDataOff)

  def _getDataPPratio(self):
    try:
      return self.res[self._dsNamePrefix+'_PumpedRatio']
    except:
      self.res[self._dsNamePrefix+'_PumpedRatio'] = \
        self.timeDelay.ones()*\
	(self.data*self.pumpOn*self.probeOn / self.dataOff.mean())
      return self.res[self._dsNamePrefix+'_PumpedRatio']
  dataPPratio = property(_getDataPPratio)
    
  def _getDataNormPPratio(self):
    try:
      return self.res[self._dsNamePrefix+'_Norm']
    except:
      self.res[self._dsNamePrefix+'_Norm'] = \
        self.timeDelay.ones()*\
	(self.data*self.pumpOn*self.probeOn/self.monitor / self.dataOffNorm.mean())
      return self.res[self._dsNamePrefix+'_Norm']
  dataNormPPratio = property(_getDataNormPPratio)

  #def _getCorrfun(self,**kwargs):
    #if self.nonlinCorrObj is None:
      #if self._nonlinCorrFile is not None:
	#cn = CorrNonLin(fina=self._nonlincorrFile)
	#corrNL = cn.getCorrFunc(**kwargs)

      #else:
	#corrNL = lambda x,im=1: x
    #else:
      #if self._corrNL is None:
        #self._corrNL = self.nonlinCorrObj.getCorrFunc()
      #corrNL = self._corrNL	
    #return corrNL
  #corrNL = property(_getCorrfun)

  #def getNonlinCorr(self,source='off'):
    #if source=='off':
      #self.nonlinCorrObj = CorrNonLin(self.dataOff)

      #data_corr = corrNL(d.cspadAzAv.data)
      #AzAvloff = loff*data_corr/d.cspadAzAv['I0']
      #iloff = AzAvloff[:,:]e
      #d.cspadAzAv['loff_corr'] = AzAvloff
      ##return iloff
      #d.cspadAzAv['MedianOff_corr'] = [tools.nanmedian(tiloff,axis=0) for tiloff in iloff]
      #d.cspadAzAv['TRDelta_corr'] = (lon*data_corr/d.cspadAzAv['I0'])-d.cspadAzAv['MedianOff_corr']
      #d.cspadAzAv['TRRatio_corr'] = (lon*data_corr/d.cspadAzAv['I0'])/d.cspadAzAv['MedianOff_corr']
      #d.cspadAzAv['norm_corr'] = (lon*data_corr/d.cspadAzAv['I0'])
      #d.cspadAzAv['TRsortedDelta_corr'] =  d.timeTool['time_binned'].ones() * d.cspadAzAv['TRDelta_corr'] 
      #d.cspadAzAv['TRsortedRatio_corr'] =  d.timeTool['time_binned'].ones() * d.cspadAzAv['TRRatio_corr'] 

    
	


    

    #if not 'I0' in d.cspadAzAv._get_keys():
      #d.cspadAzAv['I0_dat'] = nansum(d.cspadAzAv.data)
      #print "Extracting I0 from radial profile..."
      #d.cspadAzAv['I0_dat'].evaluate()
      #d.cspadAzAv['I0'] = d.cspadAzAv['I0_dat'].get_memdata()[0]
    #d.timeTool['tt_s'] = polyval(p2_tt,d.timeTool.xpp_timetool_fltpos_value)
    #ttbaseflt = d.timeTool.xpp_timetool_fltpos_value.filter([0,1000]).ones()
    #i0ttflt = ixppy.corrFilt(d.ipm2.sum,ttbaseflt*d.timeTool['tt_s'])[0]
    #i0flt = d.ipm2.sum.filter().ones()
    #gooddata = ixppy.corrFilt(d.ipm2.sum,i0ttflt*d.cspadAzAv['I0'])[0]
    #gooddata = (gooddata*d.timeTool['tt_s']).filter().ones()
    #lon = gooddata*d.eventCode.code_90.filter(True).ones()
    #tscan = d.scan['lxt_vitara_ttc']
    #d.timeTool['time'] = lon*d.timeTool['tt_s'] + tscan
    #d.timeTool['time_binned'] = d.timeTool['time'].digitize(np.arange(np.min(tscan),np.max(tscan),timeBinSize))
    #loff = gooddata*d.eventCode.code_91.filter(True).ones()
    #AzAvloff = loff*d.cspadAzAv.data/d.cspadAzAv['I0']
    #iloff = AzAvloff[:,:]
    #d.cspadAzAv['loff'] = AzAvloff
    ##return iloff
    #d.cspadAzAv['MedianOff'] = [tools.nanmedian(tiloff,axis=0) for tiloff in iloff]
    #d.cspadAzAv['TRDelta'] = (lon*d.cspadAzAv.data/d.cspadAzAv['I0'])-d.cspadAzAv['MedianOff']
    #d.cspadAzAv['TRRatio'] = (lon*d.cspadAzAv.data/d.cspadAzAv['I0'])/d.cspadAzAv['MedianOff']
    #d.cspadAzAv['norm'] = (lon*d.cspadAzAv.data/d.cspadAzAv['I0'])
    #d.cspadAzAv['TRsortedDelta'] =  d.timeTool['time_binned'].ones() * d.cspadAzAv['TRDelta'] 
    #d.cspadAzAv['TRsortedRatio'] =  d.timeTool['time_binned'].ones() * d.cspadAzAv['TRRatio'] 
#
    #d.timeTool.xpp_timetool_fltpos_value*d.ipm3.sum.filter([.3,3]).ones()*d.eventCode.code_91.filter(False).ones()

    #d.save()

    #if timeDelay is not None:
      #if delayBin is not None:
	#if type(delayBin) is np.ndarray:
	  #edges = delayBin
	#elif type(delayBin) is dict:
	  #edges = getEdges(timeDelay.R,delayBin)
	#td = timeDelay.digitize(edges)
      #else:
	#td = timeDelay
      #res = td.ones()*diffRatio
    #else:
      #res = diffRatio
    #return res
#
    #return d
