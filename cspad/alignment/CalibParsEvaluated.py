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
import CalibPars        as calp
import CalibParsDefault as cpdef
import CSPadConfigPars  as ccp

import GlobalGraphics   as gg # For test purpose in main only
import HDF5Methods      as hm # For test purpose in main only

#---------------------
#  Class definition --
#---------------------

class CalibParsEvaluated (object) :
    """This class provides access to the evaluated calibration parameters
    """
#---------------------

    def __init__ (self) :

        self.run = None # 10
        
        self.list_of_eval_types =[  'center_global'
                                   ,'rotation_index'
                                 ]
        self.mode = 'DEFAULT'

        self.setCalibParsEvaluatedFromDefault()
        self.evaluateCSPadGeometry()

#---------------------

    def setCalibParsEvaluatedFromDefault (self) :
        #print 'Set default calibration parameters for evaluated'

        self.mode = 'DEFAULT'

        self.evalpars = {}
        for type in self.list_of_eval_types :
            self.evalpars[type] = self.getCalibParsDefault (type)

#---------------------

    def getCalibParsDefault (self, type) :
        return cpdef.calibparsdefault.getCalibParsDefault (type)

#---------------------

    def setCalibParsEvaluated (self) :
        print 'Set the calibration parameters evaluated'

        self.mode = 'EVALUATED'

        # For now the only type of evaluated parameters is the 'center_global'
        self.evalpars['center_global'] = self.evaluateCSPadCenterGlobal() # np.array()
        #self.printCalibParsEvaluated ('center_global') 

        self.evaluateCSPadGeometry()

#---------------------

    def evaluateCSPadCenterGlobal(self) :
        """Evaluate the CSPad 'center_global' coordinates for all 2x1s and override the default version
        """
        offset         = calp.calibpars.getCalibPars ('offset')
        offset_corr    = calp.calibpars.getCalibPars ('offset_corr')
        marg_gap_shift = calp.calibpars.getCalibPars ('marg_gap_shift')
        quad_rotation  = calp.calibpars.getCalibPars ('quad_rotation')
        quad_tilt      = calp.calibpars.getCalibPars ('quad_tilt')

        margX    = marg_gap_shift[0][1]
        margY    = marg_gap_shift[1][1]
        gapX     = marg_gap_shift[0][2]
        gapY     = marg_gap_shift[1][2]
        shiftX   = marg_gap_shift[0][3]
        shiftY   = marg_gap_shift[1][3]

        dx = np.array([margX-gapX+shiftX,  margX-gapX-shiftX,  margX+gapX-shiftX,  margX+gapX+shiftX])
        dy = np.array([margY-gapY-shiftY,  margY+gapY-shiftY,  margY+gapY+shiftY,  margY-gapY+shiftY])
        dz = np.array([0, 0, 0, 0])

        xmin_quad = offset[0] + offset_corr[0] + dx 
        ymin_quad = offset[1] + offset_corr[1] + dy
        zmin_quad = offset[2] + offset_corr[2] + dz

        #print 'offset =\n',         offset
        #print 'offset_corr =\n',    offset_corr
        #print 'marg_gap_shift =\n', marg_gap_shift
        #print 'quad_rotation =\n',  quad_rotation        
        #print 'quad_tilt =\n',      quad_tilt    

        #print 'dx =', dx 
        #print 'dy =', dy 
        #print 'dz =', dz 

        #print 'xmin_quad =',xmin_quad  
        #print 'ymin_quad =',ymin_quad  
        #print 'zmin_quad =',zmin_quad  

        self.fill2x1CentersInQuads()

        xc_glob = np.zeros( (4,8), dtype=np.float32 )
        yc_glob = np.zeros( (4,8), dtype=np.float32 )
        zc_glob = np.zeros( (4,8), dtype=np.float32 )

        #quad_rotation_default = np.array([90,0,270,180])

        for quad in range(4) :

            #if self.mode == 'DEFAULT' : coords_in_quad = self.get2x1CentersInQuadForRotN90(quad, quad_rotation_default[quad])
            #else                      :
            coords_in_quad = self.get2x1CentersInQuadForRotN90(quad, quad_rotation[quad])

            #print 'Quad:', quad
            #print 'coords_in_quad:\n', coords_in_quad

            xc_glob[quad] = coords_in_quad[0] + xmin_quad[quad]
            yc_glob[quad] = coords_in_quad[1] + ymin_quad[quad]
            zc_glob[quad] = coords_in_quad[2] + zmin_quad[quad]

        self.center_global_evaluated = np.array([ xc_glob, yc_glob, zc_glob])

        #print 'center_global_evaluated =\n', self.center_global_evaluated
        #return self.evalpars['center_global']
        return self.center_global_evaluated

#---------------------

    def fill2x1CentersInQuads(self) :

        """fill2x1CentersInQuads(self) :"""

        marg_gap_shift = calp.calibpars.getCalibPars ('marg_gap_shift')

        quadMargX= marg_gap_shift[0][0]
        quadMargY= marg_gap_shift[1][0]
        quadMargZ= marg_gap_shift[2][0]

        center        = calp.calibpars.getCalibPars ('center')
        center_corr   = calp.calibpars.getCalibPars ('center_corr')
        quad_rotation = calp.calibpars.getCalibPars ('quad_rotation')
        quad_tilt     = calp.calibpars.getCalibPars ('quad_tilt')

        #Fill arrays [4,8] for x, y, and z:
        self.xcenter  = center[0] + center_corr[0] + quadMargX
        self.ycenter  = center[1] + center_corr[1] + quadMargY
        self.zcenter  = center[2] + center_corr[2] + quadMargZ

        self.coor_x_min = 0
        self.coor_y_min = 0
        self.coor_x_max = ccp.cspadconfig.quadDimX
        self.coor_y_max = ccp.cspadconfig.quadDimY

        #print 'self.xcenter =\n', self.xcenter  
        #print 'self.ycenter =\n', self.ycenter  
        #print 'self.zcenter =\n', self.zcenter  

#---------------------

    def get2x1CentersInQuadForRotN90(self, quad, rotN90) :
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

    def evaluateCSPadGeometry (self) :
        """Evaluate the image size and subtract xmin, ymin from all center coordinates
        """

        self.xc_arr      = self.getCalibParsEvaluated ('center_global')[0]
        self.yc_arr      = self.getCalibParsEvaluated ('center_global')[1]
        self.rot_ind_arr = self.getCalibParsEvaluated ('rotation_index')

        margin      = ccp.cspadconfig.margin
        wid2x1      = ccp.cspadconfig.wid2x1
        len2x1      = ccp.cspadconfig.len2x1 + ccp.cspadconfig.gapIn2x1
        #print 'wid2x1, len2x1 =', wid2x1, len2x1

        # 1) Find cspad x,y min and max values

        xmin,ymin,xmax,ymax=1000,1000,0,0
        for quad in range(4) :
            for segm in range(8): # loop over 0,1,2,...,7
                xc      = self.xc_arr[quad,segm]
                yc      = self.yc_arr[quad,segm]
                self.rot_ind = self.rot_ind_arr[quad,segm]
                if self.mode == 'EVALUATED' : self.rot_ind += 1 #!!!!!!!!!!!!!!! + 1 !!!!!!!!!!!!!!!!!

                #print 'rot_ind, rot_ind%2 =', rot_ind, rot_ind%2 # IT WORKS

                if self.rot_ind%2 == 0 : xsize, ysize = wid2x1/2, len2x1/2
                else                   : xsize, ysize = len2x1/2, wid2x1/2

                if xc-xsize < xmin : xmin = xc-xsize
                if xc+xsize > xmax : xmax = xc+xsize
                if yc-ysize < ymin : ymin = yc-ysize
                if yc+ysize > ymax : ymax = yc+ysize
        #print 'xmin,ymin,xmax,ymax=', xmin,ymin,xmax,ymax

        # 2) Apply offsets

        self.detDimX = int(xmax - xmin + 2*margin)
        self.detDimY = int(ymax - ymin + 2*margin)
        self.centerX = np.zeros((4, 8), dtype=np.float32)
        self.centerY = np.zeros((4, 8), dtype=np.float32)

        #print 'Image size =', self.detDimX, self.detDimY 

        for quad in range(4) :
            for segm in range(8): # loop over 0,1,2,...,7
                self.centerX[quad,segm] = self.xc_arr[quad,segm] - xmin + margin
                self.centerY[quad,segm] = self.yc_arr[quad,segm] - ymin + margin               

#---------------------

    def getCSPadGeometry (self, rotation=0) :

        rot = rotation % 4
        #print 'Returns the CSPad alignment parameters for rotation index =', rotation

        self.centerX_final   = np.zeros((4, 8), dtype=np.float32)
        self.centerY_final   = np.zeros((4, 8), dtype=np.float32)
        self.rot_index_final = np.zeros((4, 8), dtype=np.int)

        if rot % 2 == 0 :
            self.detDimX_final = self.detDimX
            self.detDimY_final = self.detDimY
        else :
            self.detDimX_final = self.detDimY
            self.detDimY_final = self.detDimX
            
        for quad in range(4) :
            for segm in range(8): # loop over 0,1,2,...,7

                self.rot_ind = self.rot_ind_arr[quad,segm]
                if self.mode == 'EVALUATED' : self.rot_ind += 1 #!!!!!!!!!!!!!!! + 1 !!!!!!!!!!!!!!!!!

                self.rot_index_final[quad,segm] = (self.rot_ind - rotation ) % 4


                if rot == 0 :
                    self.centerX_final[quad,segm] = self.centerX[quad,segm]
                    self.centerY_final[quad,segm] = self.centerY[quad,segm]

                elif rot == 1 :
                    self.centerX_final[quad,segm] =                self.centerY[quad,segm]
                    self.centerY_final[quad,segm] = self.detDimX - self.centerX[quad,segm]

                elif rot == 2 :
                    self.centerX_final[quad,segm] = self.detDimX - self.centerX[quad,segm]
                    self.centerY_final[quad,segm] = self.detDimY - self.centerY[quad,segm]

                else :
                    self.centerX_final[quad,segm] = self.detDimY - self.centerY[quad,segm]
                    self.centerY_final[quad,segm] =                self.centerX[quad,segm]

        return self.detDimX_final, self.detDimY_final, self.centerX_final, self.centerY_final, self.rot_index_final 

#---------------------

    def evaluateCSPadPixCoordinates (self, rotation=0, mirror=False) :

        dPhi        = calp.calibpars.getCalibPars ('tilt')
        detDimX, detDimY, segmX, segmY, segmRotInd = self.getCSPadGeometry (rotation-1) # !!! -1 in order to be consistent with getCSPadGeometry
        nquads      = ccp.cspadconfig.nquads
        nsects      = ccp.cspadconfig.nsects
        pixSize     = ccp.cspadconfig.pixSize

        lenCoords_um, widCoords_um = self.evaluate2x1PixCoordinates_um ()
        segmX_um = segmX*pixSize
        segmY_um = segmY*pixSize 

        detDimX_um = detDimX*pixSize
        detDimY_um = detDimY*pixSize

        print 'evaluateCSPadPixCoordinates (...):'
        print 'detDimX, detDimY =', detDimX, detDimY
        print 'detDimX_um, detDimY_um =', detDimX_um, detDimY_um
        print 'nquads, nsects   =', nquads, nsects

        self.print2DNumpyArr(segmX_um,   title='Segment center X (um):')
        self.print2DNumpyArr(segmY_um,   title='Segment center Y (um):')
        self.print2DNumpyArr(segmRotInd, title='Segment rotation index:', format='%3d')


        X, Y = np.meshgrid(lenCoords_um, widCoords_um)
        #print 'lenCoords_um.shape=', lenCoords_um.shape
        #print 'X.shape=', X.shape
        #print 'Y.shape=', Y.shape

        self.pix_global_x = np.zeros((nquads,nsects,185,388), dtype=np.float32)
        self.pix_global_y = np.zeros((nquads,nsects,185,388), dtype=np.float32)
        
        for quad in range(nquads) :
            for sect in range(nsects) :
                angle     = -(dPhi[quad][sect] + 90*segmRotInd[quad][sect] + 90)
                angle_rad = math.radians(angle)                
                S,C = math.sin(angle_rad), math.cos(angle_rad)
                Xrot, Yrot = self.rotation(X, Y, C, S)

                xc = segmX_um[quad][sect]
                yc = segmY_um[quad][sect]

                if mirror :
                    self.pix_global_x[quad][sect][:] =  Xrot + xc
                    self.pix_global_y[quad][sect][:] =  Yrot + yc
                else :
                    self.pix_global_x[quad][sect][:] = -Xrot - xc + detDimX_um-1 # mirror wrt X
                    self.pix_global_y[quad][sect][:] =  Yrot + yc
                    #self.pix_global_y[quad][sect][:] = -Yrot - yc + detDimY_um # mirror wrt Y

#---------------------

    def getCSPadPixCoordinates_um (self) :
        return self.pix_global_x, self.pix_global_y

#---------------------

    def getCSPadPixCoordinates_pix (self) :
        pixSize = ccp.cspadconfig.pixSize
        return self.pix_global_x / pixSize, self.pix_global_y / pixSize

#---------------------

    def evaluateCSPadPixCoordinatesShapedAsData(self, fname, dsname, rotation=0, mirror=False) :

        print 'Evaluate pix coordinates for fname:', fname

        self.evaluateCSPadPixCoordinates (rotation, mirror)

        ccp.cspadconfig.setCSPadConfiguration(fname, dsname, event=0)
        quadNumsInEvent  = ccp.cspadconfig.quadNumsInEvent
        indPairsInQuads  = ccp.cspadconfig.indPairsInQuads
        nquads           = ccp.cspadconfig.nquads
        nsects           = ccp.cspadconfig.nsects
        #ccp.cspadconfig.printCSPadConfigPars()

        nsects_in_data = max(indPairsInQuads.flatten()) + 1
        self.pix_global_shaped_as_data_x = np.zeros((nsects_in_data,185,388), dtype=np.float32)
        self.pix_global_shaped_as_data_y = np.zeros((nsects_in_data,185,388), dtype=np.float32)
        
        for iq in range(len(quadNumsInEvent)) :
            quad = int(quadNumsInEvent[iq]) # uint8 -> int
            for segm in range(8): # loop over ind = 0,1,2,...,7
                ind_segm_in_arr = indPairsInQuads[quad][segm]
                if ind_segm_in_arr == -1 : continue

                self.pix_global_shaped_as_data_x[ind_segm_in_arr][:] = self.pix_global_x[quad][segm][:]
                self.pix_global_shaped_as_data_y[ind_segm_in_arr][:] = self.pix_global_y[quad][segm][:]

        print 'self.pix_global_shaped_as_data_x.shape =', self.pix_global_shaped_as_data_x.shape       
        print 'evaluateCSPadPixCoordinatesShapedAsData: done'

#---------------------

    def getCSPadPixCoordinatesShapedAsData_um (self) :
        return self.pix_global_shaped_as_data_x, self.pix_global_shaped_as_data_y

#---------------------

    def getCSPadPixCoordinatesShapedAsData_pix (self) :
        pixSize = ccp.cspadconfig.pixSize
        return self.pix_global_shaped_as_data_x / pixSize, self.pix_global_shaped_as_data_y / pixSize

#---------------------

    def printCSPadPixCoordinates (self) :

        print 'self.pix_global_x =\n', self.pix_global_x
        print 'self.pix_global_y =\n', self.pix_global_y

        print 'self.pix_global_x.shape=', self.pix_global_x.shape
        print 'self.pix_global_y.shape=', self.pix_global_y.shape

#---------------------

    def rotation(self, X, Y, C, S) :
        """For numpy arryys X and Y returns the numpy arrays of Xrot and Yrot
        """
        Xrot = X*C + Y*S 
        Yrot = Y*C - X*S 
        return Xrot, Yrot

#---------------------

    def getTestImageForEntireArray(self,ds1ev) :
        """WERY SLOW METHOD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        """        
        xpix, ypix = cpeval.getCSPadPixCoordinates_pix ()
        
        #dimX,dimY = 1750, 1750
        dimX,dimY = self.detDimX, self.detDimY
        img_arr = np.zeros((dimX+1,dimY+1), dtype=np.float32)
        
        for quad in range(4) :
            for sect in range(8) :
                segm = quad*8+sect
                for row in range(185) :
                    for col in range(388) :

                        x = int( xpix[quad][sect][row][col] )
                        y = int( ypix[quad][sect][row][col] )

                        if x<0    : x=0
                        if y<0    : y=0
                        if x>dimX : x=dimX
                        if y>dimY : y=dimY

                        img_arr[x][y] = ds1ev[segm][row][col]

        return img_arr

#---------------------

    def getTestImageShapedAsData(self,ds1ev) :
        """WERY SLOW METHOD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        """        
        xpix, ypix = cpeval.getCSPadPixCoordinatesShapedAsData_pix ()

        print 'Data       array shape = ', ds1ev.shape
        print 'Coordinate array shape = ',  xpix.shape
        
        dimX,dimY = self.detDimX, self.detDimY

        img_arr = np.zeros((dimX+1,dimY+1), dtype=np.float32)

        nsect_in_arr = xpix.shape[0] 
        print 'nsect_in_arr =', nsect_in_arr

        for sect in range(nsect_in_arr) :
            for row in range(185) :
                for col in range(388) :

                    x = int( xpix[sect][row][col] )
                    y = int( ypix[sect][row][col] )

                    if x<0    : x=0
                    if y<0    : y=0
                    if x>dimX : x=dimX
                    if y>dimY : y=dimY

                    img_arr[x][y] = ds1ev[sect][row][col]

        return img_arr

#---------------------

    def evaluate2x1PixCoordinates_um (self) :

        wid2x1      = ccp.cspadconfig.wid2x1
        len2x1      = ccp.cspadconfig.len2x1
        pixSize     = ccp.cspadconfig.pixSize
        pixSizeWide = ccp.cspadconfig.pixSizeWide

        print 'wid2x1, len2x1, pixSize, pixSizeWide =', wid2x1, len2x1, pixSize, pixSizeWide

        widOffset = (184/2)*pixSize
        lenOffset = pixSizeWide + 0.5*pixSize

        lenCoords  = np.zeros((len2x1), dtype=np.float32)
        widCoords  = np.zeros((wid2x1), dtype=np.float32)

        for i in range(wid2x1) :
            widCoords[i] = -i*pixSize + widOffset

        for i in range(0,193) :
            dl = lenOffset + i*pixSize
            lenCoords[195+i] =  dl
            lenCoords[192-i] = -dl

        lenCoords[193] = -pixSizeWide/2
        lenCoords[194] =  pixSizeWide/2

        #for i in range(wid2x1) : print 'widCoords[%3d] = %10.2f' % (i,widCoords[i]) 
        #for i in range(len2x1) : print 'lenCoords[%3d] = %10.2f' % (i,lenCoords[i]) 

        return lenCoords, widCoords # in um

#---------------------

    def print2DNumpyArr( self, arr=np.zeros((4,8), dtype=np.float32), title='', format='%10.2f' ) :
        print title,
        rows,cols = arr.shape

        for row in range(rows) :
            print '\nrow ', row, ':',
            for col in range(cols) :
                print format % (arr[row][col]),                   
        print '\n'

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

        if type in self.list_of_eval_types :
            return self.evalpars[type]
        else :
            print  'WARNING: THE REQUESTED TYPE OF CALIBRATION PARS "', type, \
                   '" IS NOT FOUND IN THE AVAILABLE LIST:\n', self.list_of_eval_types
            return None

#----------------------------------------------

cpeval = CalibParsEvaluated() # Sets the default calibration parameters.

#----------------------------------------------
# In case someone decides to run this module --
#----------------------------------------------

def main_test() :

    cpeval.printCalibParsEvaluated()
    cpeval.printListOfEvaluatedTypes()
    cpeval.printCalibParsEvaluated('center_global')
    cpeval.printCalibParsEvaluated('rotation_index')

    print calp.calibpars.getCalibPars('center')
    calp.calibpars.printCalibPars() # prints the calib pars from files, if found
    print 'center_global =\n', cpeval.getCalibParsEvaluated('center_global')

if __name__ == "__main__" :

    main_test()
    sys.exit ( 'End of job' )

#----------------------------------------------
