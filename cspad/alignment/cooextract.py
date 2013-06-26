import CalibPars          as calp
import CalibParsEvaluated as cpe
import CSPadConfigPars    as ccp
import numpy as np


def getCsPadPixCoordinates(path_calib='/reg/neh/home1/dubrovin/LCLS/CSPadAlignment-v01/calib-xpp-2013-01-29', 
                           rotation=0, 
                           mirror=False):

  calp.calibpars.setCalibParsForPath (path=path_calib)
  cpe.cpeval.evaluateCSPadPixCoordinates (rotation, mirror)

  #ccp.cspadconfig.setCSPadConfiguration(fname, dsname, event=0)
  quadNumsInEvent  = np.arange(4)
  indPairsInQuads  = np.arange(32).reshape(4,-1)
  nquads           = 4
  nsects           = 8

  nsects_in_data = max(indPairsInQuads.flatten()) + 1
  x = np.zeros((nsects_in_data,185,388), dtype=np.float32)
  y = np.zeros((nsects_in_data,185,388), dtype=np.float32)
  
  for iq in range(len(quadNumsInEvent)) :
      quad = int(quadNumsInEvent[iq]) # uint8 -> int
      for segm in range(8): # loop over ind = 0,1,2,...,7
          ind_segm_in_arr = indPairsInQuads[quad][segm]
          if ind_segm_in_arr == -1 : continue

          x[ind_segm_in_arr][:] = cpe.cpeval.pix_global_x[quad][segm][:]
          y[ind_segm_in_arr][:] = cpe.cpeval.pix_global_y[quad][segm][:]
  return x,y


