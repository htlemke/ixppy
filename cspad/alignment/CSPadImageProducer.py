#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module CSPadImageProducer...
#
#------------------------------------------------------------------------

"""This module provides access to the calibration parameters

This software was developed for the SIT project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@see RelatedModule

@version $Id: 2008-09-27$

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
import scipy.ndimage as spi # rotate(...)

import CalibPars          as calp
import CalibParsEvaluated as cpe
import CSPadConfigPars    as ccp

import GlobalGraphics as gg # For test purpose in main only
import HDF5Methods    as hm # For test purpose in main only

#---------------------
#  Class definition --
#---------------------

class CSPadImageProducer (object) :
    """This class produces the device dependent CSPad image"""

#---------------------

    def getImageArrayForPair( self, arr1ev, pairNum=None ):
        """Returns the image array for pair of ASICs"""
        if pairNum == None :
            self.pair = ccp.cspadconfig.cspadPair
        else :
            self.pair = pairNum

        asic2x1 = arr1ev[self.pair,...]
        asics   = np.hsplit(asic2x1,2)
        arrgap = np.zeros( (185,3), dtype=np.int16 )
        arr2d  = np.hstack((asics[0],arrgap,asics[1]))
        return arr2d

#---------------------

    def getImageArrayForQuad( self, arr1ev, quadNum=None ):
        """Returns the image array for one quad"""

        if ccp.cspadconfig.isCSPad2x2 : # For CSpad2x2
            return self.getImageArrayForCSpad2x2Element( arr1ev )

        if quadNum == None :
            self.quad = ccp.cspadconfig.cspadQuad
        else :
            self.quad = quadNum            

        indPairsInQuads  = ccp.cspadconfig.indPairsInQuads
        pairInQaudOriInd = ccp.cspadconfig.pairInQaudOriInd

        pairXInQaud = calp.calibpars.getCalibPars ('center')[0]
        pairYInQaud = calp.calibpars.getCalibPars ('center')[1]
        dXInQaud    = calp.calibpars.getCalibPars ('center_corr')[0]
        dYInQaud    = calp.calibpars.getCalibPars ('center_corr')[1]
        offsetX     = calp.calibpars.getCalibPars ('marg_gap_shift')[0][0]
        offsetY     = calp.calibpars.getCalibPars ('marg_gap_shift')[1][0]
        dPhi        = calp.calibpars.getCalibPars ('tilt')

        #arr2dquad = np.zeros( (850,850), dtype=np.float32 ) # np.int16
        arr2dquad = np.zeros( (ccp.cspadconfig.quadDimX,ccp.cspadconfig.quadDimY), dtype=np.float32 ) # dtype=np.int16 
        #print 'arr2dquad.shape=',arr2dquad.shape

#       for ind in xrange(1): # loop over ind = 0,1,2,...,7
        for ind in range(8): # loop over ind = 0,1,2,...,7

            pair = indPairsInQuads[self.quad][ind]
            #print 'quad,ind,pair=', self.quad, ind, pair
            if pair == -1 : continue

            asic2x1 = self.getImageArrayForPair( arr1ev, pair )
            rotarr2d_0 = np.rot90(asic2x1,pairInQaudOriInd[self.quad][ind])
            #print 'rotarr2d_0.shape=',rotarr2d_0.shape
            #print 'rotarr2d.base is asic2x1 ? ',rotarr2d.base is asic2x1 

            rotarr2d = rotarr2d_0

            ixOff  = offsetX + pairXInQaud[self.quad][ind] + dXInQaud[self.quad][ind] 
            iyOff  = offsetY + pairYInQaud[self.quad][ind] + dYInQaud[self.quad][ind]
            #print 'ixOff, iyOff :', ixOff, iyOff

            # 0:185, 0:388 -> 185x391
            rot_index = pairInQaudOriInd[self.quad][ind] 

            offS = 0.5*185
            offL = 0.5*(388+3)
            #print 'offS, offL :', offS, offL

            if rot_index % 2 == 0 : self.lx0, self.ly0 = offS, offL 
            else                  : self.lx0, self.ly0 = offL, offS  

            ixOff -= self.lx0  
            iyOff -= self.ly0  

            #-------- Apply tilt angle of 2x1 sensors
            if self.tiltIsOn :
            #if ccp.confpars.cspadApplyTiltAngle :

                r0      = math.sqrt( self.lx0*self.lx0 + self.ly0*self.ly0 )
                sinPhi  = self.ly0 / r0
                cosPhi  = self.lx0 / r0

                angle  = dPhi[self.quad][ind]
                rotarr2d = spi.rotate(rotarr2d_0, angle, reshape=True, output=np.float32 )
                dimX0,dimY0 = rotarr2d_0.shape

                rdphi = r0 * abs(math.radians(angle))
                #print 'rdphi :',rdphi

                ixOff -= rdphi * sinPhi
                iyOff -= rdphi * cosPhi

                #print 'Tilt offset dx, dy=', rdphi * sinPhi, rdphi * cosPhi

            #-------- 

            ixOff = int( ixOff )
            iyOff = int( iyOff )

            dimX, dimY = rotarr2d.shape
            #print 'ixOff, iyOff =', ixOff, iyOff,           
            #print ' dimX,  dimY =', dimX, dimY           
            
            arr2dquad[ixOff:dimX+ixOff, iyOff:dimY+iyOff] += rotarr2d[0:dimX, 0:dimY]

        #print 'arr2dquad=\n', arr2dquad
        return arr2dquad

#---------------------

    def getImageArrayForCSPadElement( self, arr1ev ):
        """Returns the image array for the CSPad detector for dataset CSPadElement"""

        quadInDetOriInd = ccp.cspadconfig.quadInDetOriInd
        quadNumsInEvent = ccp.cspadconfig.quadNumsInEvent

        margX    = calp.calibpars.getCalibPars ('marg_gap_shift')[0][1]
        margY    = calp.calibpars.getCalibPars ('marg_gap_shift')[1][1]
        gapX     = calp.calibpars.getCalibPars ('marg_gap_shift')[0][2]
        gapY     = calp.calibpars.getCalibPars ('marg_gap_shift')[1][2]
        shiftX   = calp.calibpars.getCalibPars ('marg_gap_shift')[0][3]
        shiftY   = calp.calibpars.getCalibPars ('marg_gap_shift')[1][3]
        offX     = calp.calibpars.getCalibPars ('offset')[0] + calp.calibpars.getCalibPars ('offset_corr')[0] + margX
        offY     = calp.calibpars.getCalibPars ('offset')[1] + calp.calibpars.getCalibPars ('offset_corr')[1] + margY

        #self.arr2dCSpad = np.zeros( (1710,1710), dtype=np.int16 )
        #self.arr2dCSpad = np.zeros( (1750,1750), dtype=np.int16 )
        #self.arr2dCSpad = np.zeros( (1765,1765), dtype=np.float32 )
        self.arr2dCSpad = np.zeros( (ccp.cspadconfig.detDimX,ccp.cspadconfig.detDimY), dtype=np.float32 )

        quadXOffset = [offX[0]-gapX+shiftX, offX[1]-gapX-shiftX, offX[2]+gapX-shiftX, offX[3]+gapX+shiftX]
        quadYOffset = [offY[0]-gapY-shiftY, offY[1]+gapY-shiftY, offY[2]+gapY+shiftY, offY[3]-gapY+shiftY]

        #print 'quadXOffset = ', quadXOffset 
        #print 'quadYOffset = ', quadYOffset 

        #for iq in range(1) :
        for iq in range(len(quadNumsInEvent)) :
            quad = int(quadNumsInEvent[iq]) # uint8 -> int
            arr2dquad = self.getImageArrayForQuad(arr1ev, quad)
            rotarr2d = np.rot90(arr2dquad,quadInDetOriInd[quad])
            #print 'rotarr2d.shape=',rotarr2d.shape
            dimX,dimY = rotarr2d.shape

            ixOff = quadXOffset[quad]
            iyOff = quadYOffset[quad]

            self.arr2dCSpad[ixOff:dimX+ixOff, iyOff:dimY+iyOff] += rotarr2d[0:dimX, 0:dimY]

        return self.arr2dCSpad

#---------------------

    def getImageArrayForDet( self, arr1ev ):
        """Returns the image array for entire CSpad detector"""       

        if cp.confpars.eventCurrent == self.eventWithAlreadyGeneratedCSpadDetImage :
            #print 'Use already generated image for CSpad and save time'
            return self.arr2dCSpad

        if ccp.cspadconfig.isCSPad2x2 : # For CSpad2x2
            self.arr2dCSpad = self.getImageArrayForCSpad2x2Element( arr1ev )

        else : # For regular CSPad detector
            self.arr2dCSpad = self.getImageArrayForCSPadElement( arr1ev )

        self.eventWithAlreadyGeneratedCSpadDetImage = cp.confpars.eventCurrent

        if cp.confpars.bkgdSubtractionIsOn : self.arr2dCSpad -= cp.confpars.arr_bkgd
        if cp.confpars.gainCorrectionIsOn  : self.arr2dCSpad *= cp.confpars.arr_gain

        return self.arr2dCSpad

#---------------------

    def getImageArrayForCSpad2x2ElementPair( self, arr1ev, pairNum=None ):
        """Returns the image array for pair of ASICs"""
        if pairNum == None :
            self.pair = ccp.cspadconfig.cspadPair
        else :
            self.pair = pairNum

        #arr2x1 = arr1ev[0:185,0:388,self.pair]
        arr2x1 = arr1ev[:,:,self.pair]
        asics  = np.hsplit(arr2x1,2)
        arrgap = np.zeros ((185,3), dtype=np.float32)
        arr2d  = np.hstack((asics[0],arrgap,asics[1]))
        return arr2d

#---------------------

    def getImageArrayForCSpad2x2Element( self, arr1ev ):
        """Returns the image array for the CSpad2x2Element or CSpad2x2"""       

        arr2x1Pair0 = self.getImageArrayForCSpad2x2ElementPair(arr1ev,0)
        arr2x1Pair1 = self.getImageArrayForCSpad2x2ElementPair(arr1ev,1)
        wid2x1      = arr2x1Pair0.shape[0]
        len2x1      = arr2x1Pair0.shape[1]

        arrgapV = np.zeros( (20,len2x1), dtype=np.float ) # dtype=np.int16 
        arr2d   = np.vstack((arr2x1Pair0, arrgapV, arr2x1Pair1))

        #print 'arr2d.shape=', arr2d.shape
        #print 'arr2d=',       arr2d
        return arr2d

#---------------------
#---------------------
#---------------------
# New approach to the CSPad geometry
# All 2x1 centers of entire CSPad are defined in the same coordinate frame.
#---------------------
#---------------------
#---------------------

    def getCSPadArrayWithGap(self, arr, gap=3) :
        #print 'getCSPadArrayWithGap(...): Input array shape =', arr.shape
        arr.shape = (arr.size/388,388)
        nrows,ncols = arr.shape # (32*185,388) = (5920,388) # <== expected input array shape for all sections
        if ncols != 388 or nrows<185 :
            print 'getCSPadArrayWithGap(...): WARNING! UNEXPECTED INPUT ARRAY SHAPE =', arr.shape
            return arr
        arr_gap = np.zeros( (nrows,gap), dtype=np.int16 )
        arr_halfs = np.hsplit(arr,2)
        arr_with_gap = np.hstack((arr_halfs[0], arr_gap, arr_halfs[1]))
        arr_with_gap.shape = (nrows,ncols+gap)
        return arr_with_gap

#---------------------

    def getCSPadImage(self, arr) :

        #pairInQaudOriInd = ccp.cspadconfig.pairInQaudOriInd
        #quadInDetOriInd  = ccp.cspadconfig.quadInDetOriInd
        quadNumsInEvent  = ccp.cspadconfig.quadNumsInEvent
        indPairsInQuads  = ccp.cspadconfig.indPairsInQuads
        gap              = ccp.cspadconfig.gapIn2x1
        dPhi             = calp.calibpars.getCalibPars ('tilt')

        #print 'quadNumsInEvent =', quadNumsInEvent
        #print 'indPairsInQuads =\n',  indPairsInQuads

        arr_all       = self.getCSPadArrayWithGap(arr, gap)

        dim3 = arr_all.shape[arr.ndim-1] # should be 388 or 388+gap
        dim2 = 185
        dim1 = arr_all.size/dim3/dim2
        arr_all.shape = (dim1,dim2,dim3) # Reshape for quad and segment indexes

        #print 'dims          =', dim1, dim2, dim3
        #print 'arr_all.shape =', arr_all.shape

        arr_cspad_img = np.zeros( (self.detDimX,self.detDimY), dtype=np.float32 )

        for indq in range(len(quadNumsInEvent)) :
            quad = int(quadNumsInEvent[indq]) # uint8 -> int

            for segm in range(8): # loop over 0,1,2,...,7

                ind_segm_in_arr = indPairsInQuads[quad][segm]
                #print 'quad, segm, ind_segm_in_arr =', quad, segm, ind_segm_in_arr
                if ind_segm_in_arr == -1 : continue

                #arr_segm = arr_all[quad,segm,:] # old slyle of shape...
                arr_segm = arr_all[ind_segm_in_arr,:]    

                if self.tiltIsOn :
                    angle  = dPhi[quad][segm] + 90*self.segmRotInd[quad][segm]
                    #print 'angle = %f ' % (angle) 
                    #arr_segm_rot0 = np.rot90(arr_segm, self.segmRotInd[quad][segm])
                    arr_segm_rot  = spi.rotate(arr_segm, angle, reshape=True, output=np.float32 )
                else :
                    arr_segm_rot = np.rot90(arr_segm, self.segmRotInd[quad][segm])
        
                nrows, ncols = arr_segm_rot.shape
                #print 'quad, segm, nrows, ncols = ', quad, segm, nrows, ncols
        
                xOff = self.segmX[quad][segm] - nrows/2
                yOff = self.segmY[quad][segm] - ncols/2
        
                arr_cspad_img[xOff:nrows+xOff, yOff:ncols+yOff] += arr_segm_rot[0:nrows, 0:ncols]

        if self.mirror : return  np.fliplr(arr_cspad_img)
        else           : return  arr_cspad_img

#---------------------

    def printInputPars(self) :
        print '\nCSPadImageProducer(): printInputPars()'
        print 'self.rotation =', self.rotation
        print 'self.mirror   =', self.mirror
        print 'self.tiltIsOn =', self.tiltIsOn

#---------------------

    def printGeometryPars(self) :
        print '\nCSPadImageProducer(): printGeometryPars()'
        print 'self.detDimX, self.detDimY =', self.detDimX, self.detDimY
        print 'segmX =\n',      self.segmX
        print 'segmY =\n',      self.segmY
        print 'segmRotInd =\n', self.segmRotInd

#---------------------

    def __init__ (self, rotation=0, tiltIsOn=False, mirror=False) :
        #print 'CSPadImageProducer(): Initialization'

        self.rotation = rotation
        self.mirror   = mirror
        self.tiltIsOn = tiltIsOn

        self.detDimX, self.detDimY, self.segmX, self.segmY, self.segmRotInd = \
                      cpe.cpeval.getCSPadGeometry (rotation)  

        #self.printInputPars()
        #self.printGeometryPars()

#----------------------------------------------

def main_calib() :

    print 'Start test in main_calib()'

    #calp.calibpars.setCalibPars( run      = 9,
    #                             calibdir = '/reg/d/psdm/CXI/cxi35711/calib',
    #                             group    = 'CsPad::CalibV1',
    #                             source   = 'CxiDs1.0:Cspad.0' )

    runnum=0

    #path_calib = '/reg/neh/home/dubrovin/LCLS/CSPadAlignment-v01/calib-cxi37411-r0039-Dsd/' 
    #path_calib = '/reg/neh/home/dubrovin/LCLS/CSPadAlignment-v01/calib-cxi37411-r0080-Ds1' 
    #path_calib = '/reg/neh/home/dubrovin/LCLS/CSPadAlignment-v01/calib-cxi35711-r0009-det' 
    #path_calib = '/reg/d/psdm/CXI/cxi37411/calib/CsPad::CalibV1/CxiDs1.0:Cspad.0' 
    #path_calib = '/reg/d/psdm/CXI/cxi35711/calib/CsPad::CalibV1/CxiDs1.0:Cspad.0' 
    #path_calib = '/reg/d/psdm/CXI/cxi37411/calib/CsPad::CalibV1/CxiDsd.0:Cspad.0' 
    path_calib = '/reg/d/psdm/xpp/xpp47712/calib/CsPad::CalibV1/XppGon.0:Cspad.0' 
    
    #fname  = '/reg/d/psdm/CXI/cxi35711/hdf5/cxi35711-r0009.h5'
    #fname  = '/reg/d/psdm/CXI/cxi37411/hdf5/cxi37411-r0080.h5'
    #fname  = '/reg/d/psdm/CXI/cxi37411/hdf5/cxi37411-r0039.h5'
    #fname  = '/reg/d/psdm/XPP/xpp47712/hdf5/xpp47712-r0043.h5'
    #fname  = '/reg/d/psdm/CXI/cxi80410/hdf5/cxi80410-r0628.h5'
    fname, runnum = '/reg/d/psdm/XPP/xppcom10/hdf5/xppcom10-r1437.h5', 1437

    #dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/CxiDs1.0:Cspad.0/data'
    #dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/CxiDsd.0:Cspad.0/data'
    dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/XppGon.0:Cspad.0/data'

    event  = 0

    print 'Load calibration parameters from', path_calib 
    calp.calibpars.setCalibParsForPath ( run=runnum, path=path_calib )

    print 'Get raw CSPad event %d from file %s \ndataset %s' % (event, fname, dsname)
    ds1ev = hm.getOneCSPadEventForTest( fname, dsname, event )
    print 'ds1ev.shape = ',ds1ev.shape

    print 'Make the CSPad image from raw array'
    #cspadimg = CSPadImageProducer(rotation=3, tiltIsOn=False, mirror=False)
    cspadimg = CSPadImageProducer(rotation=0, tiltIsOn=False, mirror=False)
    cspadimg.printInputPars()
    cspadimg.printGeometryPars()
    #arr = cspadimg.getImageArrayForPair( ds1ev, pairNum=3 )
    #arr = cspadimg.getImageArrayForQuad( ds1ev, quadNum=2 )
    #arr = cspadimg.getImageArrayForCSPadElement( ds1ev )
    arr = cspadimg.getCSPadImage( ds1ev )
    #print 'arr = \n',arr

    AmpRange = (1700,2000)
    #AmpRange = (0, 100)

    print 'Plot CSPad image'
    gg.plotImage(arr,range=AmpRange,figsize=(11.6,10))
    gg.move(200,100)
    #gg.plotImageAndSpectrum(arr,range=(1,2001))
    gg.plotSpectrum(arr,range=AmpRange)
    gg.move(50,50)
    print 'To EXIT the test click on "x" in the top-right corner of each plot window.'
    gg.show()

#----------------------------------------------

if __name__ == "__main__" :

    main_calib()
    sys.exit ( 'End of test.' )

#----------------------------------------------
