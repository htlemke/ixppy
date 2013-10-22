#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module CalibParsEvaluated...
#
#------------------------------------------------------------------------

"""This module provides access to the calibration parameters

This software was developed for the SIT project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@see RelatedModule

@version $Id: 2008-09-22$

@author Mikhail S. Dubrovin
"""

#------------------------------
#  Module's version from CVS --
#------------------------------
__version__ = "$Revision: 4 $"
# $Source$

#----------
#  Imports 
#----------
import sys
import os
import math
import numpy as np
import CalibParsDefault as cpd
import CSPadConfigPars  as ccp

#---------------------
#  Class definition --
#---------------------

class CSPADCalibParsEvaluated (object) :
    """This class evaluates calibration parameters using other types of parameters

       Currently it evaluates the center_global calibration parameters only.
       This class should not be used directly in any application. Use CalibPars() for all needs.
       If necessary, it will be instatiated and used automatically in CalibPars().

       Interface
       =========

       import CalibPars as calp

       path='/reg/d/psdm/xpp/xpptut13/calib/CsPad::CalibV1/XppGon.0:Cspad.0/'
       calibpars = calp.CalibPars(path, run=123)
       cpeval = CSPADCalibParsEvaluated (calibpars)
       center_global = cpeval.getCalibParsEvaluated('center_global') # np.array of shape=(3,8,4)

       Test methods
       ------------
       cpeval.printListOfEvaluatedTypes()
       cpeval.printCalibParsEvaluated() # for all evaluated parameters
       cpeval.printCalibParsEvaluated('center_global')
    """

    list_of_eval_types =['center_global']

#---------------------

    def __init__ (self, calibpars=None) :

        self.calibpars = calibpars

        self.evalpars = {}
        self.evalpars_status = {}

        self.setCalibParsDefault()
        #self.setCalibParsEvaluated()

#---------------------

    def setCalibParsDefault (self) :
        """Sets default calibration parameters for all types included in self.list_of_eval_types
        """
        for type in self.list_of_eval_types :
            self.evalpars       [type] = self.getCalibParsDefault (type)
            self.evalpars_status[type] = 'DEFAULT'

#---------------------

    def getCalibParsDefault (self, type) :
        """Local call returning default parameters.
           Default parameters are available through the singleton object cpd.calibparsdefault.
        """
        return cpd.calibparsdefault.getCalibParsDefault (type)

#---------------------

    def setCalibParsEvaluated (self) :
        """Sets evaluated or default parameters at instatiation.
        """
        try : 
            self.evalpars       ['center_global'] = self.evaluateCSPadCenterGlobal() # np.array()
            self.evalpars_status['center_global'] = 'EVALUATED'
        except :
            self.getCalibParsDefault ('center_global')

#---------------------

    def evaluateCSPadCenterGlobal(self) :
        """Evaluate the CSPAD 'center_global' coordinates for all 2x1s.

        1. 2x1 center coordinates are evaluated in the MATRIX coordinate system,
        with origin in the top left corner. Each quad should be properly rotated.
        2. All 2x1 center coordinates are rotated to the OPICAL MEASUREMENT coordinate system.
        """
        offset         = self.calibpars.getCalibPars ('offset')
        offset_corr    = self.calibpars.getCalibPars ('offset_corr')
        marg_gap_shift = self.calibpars.getCalibPars ('marg_gap_shift')
        quad_rotation  = self.calibpars.getCalibPars ('quad_rotation')
        quad_tilt      = self.calibpars.getCalibPars ('quad_tilt')

        #quadMargX, margX, gapX, shiftX = marg_gap_shift[0,:]
        #quadMargY, margY, gapY, shiftY = marg_gap_shift[1,:]
        #quadMargZ, margZ, gapZ, shiftZ = marg_gap_shift[2,:]

        margX,  margY,  margZ  = marg_gap_shift[:,1]
        gapX,   gapY,   gapZ   = marg_gap_shift[:,2]
        shiftX, shiftY, shiftZ = marg_gap_shift[:,3]

        dx = np.array([margX-gapX+shiftX,  margX-gapX-shiftX,  margX+gapX-shiftX,  margX+gapX+shiftX])
        dy = np.array([margY-gapY-shiftY,  margY+gapY-shiftY,  margY+gapY+shiftY,  margY-gapY+shiftY])
        dz = np.array([0, 0, 0, 0])

        xmin_quad = offset[0] + offset_corr[0] + dx 
        ymin_quad = offset[1] + offset_corr[1] + dy
        zmin_quad = offset[2] + offset_corr[2] + dz

        self.fill2x1CentersInQuads()

        xc_glob = np.zeros( (4,8), dtype=np.float32 )
        yc_glob = np.zeros( (4,8), dtype=np.float32 )
        zc_glob = np.zeros( (4,8), dtype=np.float32 )

        quad_rotation = np.array([180, 90, 0, 270]) # Rotation of quads in MATRIX coordinate system
        #quad_rotation = np.array([90, 0, 270, 180]) # Rotation of quads in OPTICAL coordinate system
        #print 'quad_rotation', quad_rotation
        #print 'offset\n', offset

        for quad in range(4) :

            coords_in_quad = self.get2x1CentersInQuadForRotN90(quad, quad_rotation[quad])

            xc_glob[quad] = coords_in_quad[0] + xmin_quad[quad]
            yc_glob[quad] = coords_in_quad[1] + ymin_quad[quad]
            zc_glob[quad] = coords_in_quad[2] + zmin_quad[quad]

        xc_glob, yc_glob, zc_glob = self.get2x1CentersInDetForRot270(xc_glob, yc_glob, zc_glob) # Transformation from MATRIX to OPTICAL coordinate system

        self.center_global_evaluated = np.array([ xc_glob, yc_glob, zc_glob])

        #print 'center_global_evaluated =\n', self.center_global_evaluated
        #return self.evalpars['center_global']
        return self.center_global_evaluated

#---------------------

    def fill2x1CentersInQuads(self) :
        """Evaluates 2x1 center coordinates in quads.
        """

        marg_gap_shift = self.calibpars.getCalibPars ('marg_gap_shift')
        center         = self.calibpars.getCalibPars ('center')
        center_corr    = self.calibpars.getCalibPars ('center_corr')

        quadMargX,  quadMargY,  quadMargZ  = marg_gap_shift[:,0]

        #Fill arrays [4,8] for x, y, and z:
        self.xcenter  = center[0] + center_corr[0] + quadMargX
        self.ycenter  = center[1] + center_corr[1] + quadMargY
        self.zcenter  = center[2] + center_corr[2] + quadMargZ

        self.coor_x_min = 0
        self.coor_y_min = 0
        self.coor_x_max = ccp.cspadconfig.quadDimX
        self.coor_y_max = ccp.cspadconfig.quadDimY

        #print 'center\n', center
        #print 'self.xcenter =\n', self.xcenter  
        #print 'self.ycenter =\n', self.ycenter  
        #print 'self.zcenter =\n', self.zcenter  

#---------------------

    def get2x1CentersInDetForRot270(self, xc_glob, yc_glob, zc_glob) :
        """Rotation of the 2x1 center coordinates from MATRIX to OPTICAL coordinate system.
        """

        self.det_x_min = 0
        self.det_y_min = 0
        self.det_x_max = ccp.cspadconfig.detDimX
        self.det_y_max = ccp.cspadconfig.detDimY

        xc_rot =  yc_glob - self.det_y_min
        yc_rot = -xc_glob + self.det_x_max
        zc_rot =  zc_glob

        return xc_rot, yc_rot, zc_rot

#---------------------

    def get2x1CentersInQuadForRotN90(self, quad, rotN90) :
        """Rotation method for quads.
        """
        if   rotN90 ==   0 : return self.get2x1CentersInQuadForRot000(quad)
        elif rotN90 ==  90 : return self.get2x1CentersInQuadForRot090(quad)
        elif rotN90 == 180 : return self.get2x1CentersInQuadForRot180(quad)
        elif rotN90 == 270 : return self.get2x1CentersInQuadForRot270(quad)
        #else :
        #    print 'ERROR in get2x1CentersInQuadForRotN90(self, quad, rotN90):\n',\
        #          'Invalid rotation angle=', rotN90, 'for quad=', quad
        #    sys.exit('Exit application due to error.')

#---------------------

    def get2x1CentersInQuadForRot000(self, quad) :
        return  np.array([ self.xcenter[quad] - self.coor_x_min, \
                           self.ycenter[quad] - self.coor_y_min, \
                           self.zcenter[quad] ])

#---------------------

    def get2x1CentersInQuadForRot090(self, quad) :
        return  np.array([-self.ycenter[quad] + self.coor_y_max, \
                           self.xcenter[quad] - self.coor_x_min, \
                           self.zcenter[quad] ]) 

#---------------------

    def get2x1CentersInQuadForRot180(self, quad) :
        return  np.array([-self.xcenter[quad] + self.coor_x_max, \
                          -self.ycenter[quad] + self.coor_y_max, \
                           self.zcenter[quad] ]) 

#---------------------

    def get2x1CentersInQuadForRot270(self, quad) :
        return  np.array([ self.ycenter[quad] - self.coor_y_min, \
                          -self.xcenter[quad] + self.coor_x_max, \
                           self.zcenter[quad] ])  
    
#---------------------

#    def rotation(self, X, Y, C, S) :
#        """For numpy arryys X and Y returns the numpy arrays of Xrot and Yrot
#        """
#        Xrot = X*C + Y*S 
#        Yrot = Y*C - X*S 
#        return Xrot, Yrot

#---------------------

    def printListOfEvaluatedTypes (self) :
        print 'printListOfEvaluatedTypes(): list_of_eval_types:', self.list_of_eval_types

#---------------------

    def printCalibParsEvaluatedAll (self) :
        for type in self.list_of_eval_types :
            print '\nEvaluated calibration parameter type "' + type + '" with shape', self.evalpars[type].shape
            print self.evalpars[type]

#---------------------

    def printCalibParsEvaluated (self, partype=None) :
        """Print the calibration prarameters of specified partype or all for dafault.
        """        
        if partype==None :
            for type in self.list_of_eval_types :
                print '\nprintCalibParsEvaluated(): Evaluated calibration parameter type "' + type + '" with shape', self.evalpars[type].shape
                print self.evalpars[type]
        else :
            if partype in self.list_of_eval_types :
                print '\nprintCalibParsEvaluated(): Evaluated calibration parameter type "' + partype + '" with shape', self.evalpars[partype].shape
                print self.evalpars[partype]
            else :
                print  'WARNING: THE REQUESTED TYPE OF CALIBRATION PARS "', partype, \
                       '" IS NOT FOUND IN THE AVAILABLE LIST:\n', self.list_of_eval_types
            
#---------------------

    def getCalibParsEvaluated (self, type) :
        """Evaluates the parameters for current state of the self.calibpars,
           and returns parameters for current type.
        """

        self.setCalibParsEvaluated()

        if type in self.list_of_eval_types :
            return self.evalpars[type]
        else :
            print  'WARNING: THE REQUESTED TYPE OF CALIBRATION PARS "', type, \
                   '" IS NOT FOUND IN THE AVAILABLE LIST:\n', self.list_of_eval_types
            return None

#----------------------------------------------

import CalibPars as calp # for test purpose only

def main_test() :
    """Test of the basic interface.
    """

    #path='/reg/d/psdm/CXI/cxi35711/calib/CsPad::CalibV1/CxiDs1.0:Cspad.0'
    calibpars = calp.CalibPars(path='/reg/d/psdm/xpp/xpptut13/calib/CsPad::CalibV1/XppGon.0:Cspad.0/', run=10)
    cpeval = CSPADCalibParsEvaluated (calibpars)

    cpeval.printListOfEvaluatedTypes()
    print 'center_global =\n', cpeval.getCalibParsEvaluated('center_global')

    #cpeval.printCalibParsEvaluated() # for all evaluated parameters
    #cpeval.printCalibParsEvaluated('center_global')

#---------------------

if __name__ == "__main__" :

    main_test()
    sys.exit ( 'End of job' )

#----------------------------------------------
