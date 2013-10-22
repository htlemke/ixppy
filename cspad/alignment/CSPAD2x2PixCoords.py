#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module CSPAD2x2PixCoords...
#
#------------------------------------------------------------------------

"""
This software was developed for the SIT project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@see RelatedModule

@version $Id: 2013-02-01$

@author Mikhail S. Dubrovin
"""

#--------------------------------
#  Module's version from CVS --
#--------------------------------
__version__ = "$Revision: 4 $"
# $Source$
#--------------------------------

import sys
import math
import numpy as np
from time import time

#import matplotlib.pyplot as plt

from PixCoords2x1 import *
import GlobalGraphics as gg # For test purpose in main only

#------------------------------

class CSPAD2x2PixCoords (PixCoords2x1) :
    """Class for generation of CSPad2x2 pixel coordinate array with and without data base

       Interface
       =========
       1.1 Instantiation with external geometry parameters:
 
           All parameters optional. Default values will be used if parameters are not specified.
           xc    = np.array([198., 198.]) * PixCoords2x1.pixs # 109.92 
           yc    = np.array([ 95., 308.]) * PixCoords2x1.pixs # 109.92
           tilt  = np.array([  0.,   0.])
           coord = CSPAD2x2PixCoords(xc_um=xc, yc_um=yc, tilt_deg=tilt)


       1.2 Instantiation with regular calibration parameters:

           path  = '/reg/d/psdm/mec/mec73313/calib/CsPad2x2::CalibV1/MecTargetChamber.0:Cspad2x2.1/'
           run   = 123
           calib = CSPAD2x2CalibPars(path, run)
           coord = CSPAD2x2PixCoords(calib)

       1.2 Access methods:
           Get arrays of pixel coordinates in mu with shape: [2,185,388]
           X,Y = coord.get_cspad2x2_pix_coordinate_arrays_mu ()
           or in pixels:
           X,Y = coord.get_cspad2x2_pix_coordinate_arrays_pix ()

       1.3 Get CSPAD2X2 image  
           data.shape = (185,388,2) <- preferable shape, converted if necessary...
           img = coord.get_cspad2x2_image(data)           

       2.  Global methods
       2.1 Conversion between shapes of arr2x2.shape=(185,388,2) and arrTwo2x1.shape=(2,185,388)
           arrTwo2x1 = data2x2ToTwo2x1(arr2x2)
           arr2x2    = two2x1ToData2x2(arrTwo2x1)
     """

##------------------------------

    sects = 2 # Total number of sections in quad
    xc_um_def    = np.array([198., 198.]) * PixCoords2x1.pixs # 109.92
    yc_um_def    = np.array([ 95., 308.]) * PixCoords2x1.pixs # 109.92
    tilt_deg_def = np.array([  0.,   0.])

##------------------------------

    def __init__ (sp, calib=None, xc_um=None, yc_um=None, tilt_deg=None, use_wide_pix_center=False) :
        #print 'CSPAD2x2PixCoords.__init__(...)'

        PixCoords2x1.__init__ (sp, use_wide_pix_center)

        if calib == None :

            if xc_um == None : sp.xc = sp.xc_um_def
            else             : sp.xc = xc_um

            if yc_um == None : sp.yc = sp.yc_um_def
            else             : sp.yc = yc_um

            sp.zc = 0

            sp.tilt       = tilt_deg

        else :
            [sp.xc,sp.yc,sp.zc] = calib.getCalibPars('center') * PixCoords2x1.pixs # 109.92
            sp.tilt       = calib.getCalibPars('tilt')
            #print 'USE xc, yc, zc =', sp.xc, sp.yc, sp.zc

        sp.calib = calib
        sp.make_cspad2x2_pix_coordinate_arrays (sp.xc, sp.yc, sp.tilt)

#------------------------------

    def print_cspad2x2_geometry_pars (sp) :
        print 'print_cspad2x2_geometry_pars(): xc, yc, zc, tilt =', sp.xc, sp.yc, sp.zc, sp.tilt

#------------------------------

    def make_cspad2x2_pix_coordinate_arrays (sp, xc_um, yc_um, tilt_deg=None) : # All lists of size[2]
        """Makes [2,185,388] cspad pixel x and y coordinate arrays"""        
        #sp.make_maps_of_2x1_pix_coordinates()

        sp.x_pix_um = np.zeros((sp.sects,sp.rows,sp.cols), dtype=np.float32)
        sp.y_pix_um = np.zeros((sp.sects,sp.rows,sp.cols), dtype=np.float32)

        angle_deg = [180,180]
        if tilt_deg != None : angle_deg += tilt_deg
 
        for sect in range(sp.sects) :

            angle_rad = math.radians(angle_deg[sect])                
            S,C = math.sin(angle_rad), math.cos(angle_rad)
            Xrot, Yrot = rotation(sp.x_map2x1_um, sp.y_map2x1_um, C, S)

            sp.x_pix_um[sect][:] =  Xrot + xc_um[sect]
            sp.y_pix_um[sect][:] =  Yrot + yc_um[sect]

        sp.x_pix_um -= sp.x_pix_um.min()
        sp.y_pix_um -= sp.y_pix_um.min() 

        sp.x_pix_pix = (sp.x_pix_um/sp.pixs+0.25).astype(int) 
        sp.y_pix_pix = (sp.y_pix_um/sp.pixs+0.25).astype(int)

        sp.x_pix_shapeed_as_data_pix = two2x1ToData2x2(sp.x_pix_pix)
        sp.y_pix_shapeed_as_data_pix = two2x1ToData2x2(sp.y_pix_pix)

#------------------------------

    def get_cspad2x2_pix_coordinate_arrays_um (sp) : 
        return sp.x_pix_um, sp.y_pix_um

    def get_cspad2x2_pix_coordinate_arrays_pix (sp) : 
        return sp.x_pix_pix, sp.y_pix_pix

    def get_cspad2x2_pix_coordinate_arrays_shapeed_as_data_um (sp) : 
        return two2x1ToData2x2(sp.x_pix_um), \
               two2x1ToData2x2(sp.y_pix_um)

    def get_cspad2x2_pix_coordinate_arrays_shapeed_as_data_pix (sp) : 
        return sp.x_pix_shapeed_as_data_pix, \
               sp.y_pix_shapeed_as_data_pix

#------------------------------

    def print_cspad2x2_coordinate_arrays(sp) :
        print 'print_cspad2x2_coordinate_arrays()'        
        print 'sp.x_pix_um:\n',      sp.x_pix_um
        print 'sp.x_pix_um.shape =', sp.x_pix_um.shape
        print 'sp.y_pix_um:\n',      sp.y_pix_um
        print 'sp.y_pix_um.shape =', sp.y_pix_um.shape

#------------------------------

    def get_cspad2x2_image(sp, data_arr=None) : # preferable data_arr.shape=(185,388,2) like in data
        """ Test of coordinate arrays, plot image map.
        """
        iX,iY = sp.get_cspad2x2_pix_coordinate_arrays_shapeed_as_data_pix ()
        data = data_arr
        if data_arr != None and data_arr.shape == (2,185,388) :
            data = two2x1ToData2x2(data_arr)
        #if data_arr != None and data_arr.shape == (185,388,2) :
        #    data = data2x2ToTwo2x1(data_arr)
        #print 'data.shape for image =', data.shape
        return gg.getImageFromIndexArrays(iX,iY,data) # All arrays should have the same shape

#------------------------------
#------------------------------
#------ Global methods --------
#------------------------------
#------------------------------

def data2x2ToTwo2x1(arr2x2) :
    """Converts array shaped as CSPAD2x2 data (185,388,2)
    to two 2x1 arrays with shape=(2,185,388)
    """
    return np.array([arr2x2[:,:,0], arr2x2[:,:,1]])

#------------------------------

def two2x1ToData2x2(arrTwo2x1) :
    """Converts array shaped as two 2x1 arrays (2,185,388)
    to CSPAD2x2 data shape=(185,388,2)
    """
    arr2x2 = np.array(zip(arrTwo2x1[0].flatten(), arrTwo2x1[1].flatten()))
    arr2x2.shape = (185,388,2)
    return arr2x2

#------------------------------
#------------------------------
#------------------------------
#----------- TEST -------------
#------------------------------
#------------------------------
#------------------------------

def test_reshaping_arrs_for_cspad2x2() :

    raw_arr = np.arange(185*388*2)
    raw_arr.shape = (185,388,2)

    ord_arr = data2x2ToTwo2x1(raw_arr)
    tst_arr = two2x1ToData2x2(ord_arr)

    print 'raw_arr:', raw_arr
    print 'ord_arr:', ord_arr
    print 'tst_arr:', tst_arr

    print 'raw_arr.shape:', raw_arr.shape
    print 'ord_arr.shape:', ord_arr.shape
    print 'tst_arr.shape:', tst_arr.shape

    if np.equal(tst_arr,raw_arr).all() : print 'Arrays are equal after two transformations'
    else                               : print 'Arrays are NOT equal after two transformations'

#------------------------------

def test_of_coord_arrs(coord, calib=None) :
    """ Test of coordinate arrays, plot image map.
    """

    #coord.print_cspad2x2_coordinate_arrays()
    #iX,iY = coord.get_cspad2x2_pix_coordinate_arrays_pix ()
    iX,iY = coord.get_cspad2x2_pix_coordinate_arrays_shapeed_as_data_pix ()

    W         = None
    amp_range = (-1, 2)
    if calib != None :
        W = calib.getCalibPars('pedestals')
        print ' W.shape =', W.shape # W.shape = (185, 388, 2 )
        amp_range = (200,600)

    print 'iX.shape =', iX.shape
    print 'iY.shape =', iY.shape

    t0_sec = time()
    #img2d = gg.getImageAs2DHist(X,Y,W=None)
    img2d = gg.getImageFromIndexArrays(iX,iY,W)
    print 'Consumed time to create image (sec) =', time()-t0_sec

    
    #gg.plotImageLarge(img2d, amp_range=(-1, 2), figsize=(12,11))
    gg.plotImageLarge(img2d, amp_range = amp_range, figsize=(12,11))
    gg.show()

#------------------------------
import PyCSPadImage.HDF5Methods as hm # For test purpose in main only

def test_of_image(coord, calib=None) :
    """ Test of coordinate arrays, plot image map.
    """
    #fname = '/reg/d/psdm/xpp/xpptut13/hdf5/xppi0513-r0008.h5'
    fname = '/reg/d/psdm/xpp/xppi0513/hdf5/xppi0513-r0008.h5'
    dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad2x2::ElementV1/XppGon.0:Cspad2x2.1/data'
    run   = 123
    dset = hm.getDataSetForOneEvent( fname, dsname, event  = 0 ) 
    iX,iY = coord.get_cspad2x2_pix_coordinate_arrays_shapeed_as_data_pix ()

    #dset = calib.getCalibPars('pedestals')
    print ' dset.shape =', dset.shape # dset.shape = (185, 388, 2 )
    t0_sec = time()
    img2d = gg.getImageFromIndexArrays(iX,iY,dset)
    print 'Consumed time to create image (sec) =', time()-t0_sec

    gg.plotImageLarge(img2d, amp_range=None, figsize=(12,11))
    gg.show()

#------------------------------

def test_instantiation_1 () :
    """ Instantiation with external geometry parameters.
    """
    xc    = np.array([198., 202.]) * PixCoords2x1.pixs # 109.92 
    yc    = np.array([ 95., 308.]) * PixCoords2x1.pixs # 109.92
    tilt  = np.array([  0.,   2.])

    t0_sec = time()
    coord = CSPAD2x2PixCoords(xc_um=xc, yc_um=yc, tilt_deg=tilt)
    coord.print_cspad2x2_geometry_pars()
    print 'Consumed time for coordinate arrays (sec) =', time()-t0_sec
    return coord

#------------------------------

from CSPAD2x2CalibPars import *

def test_instantiation_2() :
    """ Instantiation with regular calibration parameters.
    """
    #path  = '/reg/d/psdm/mec/mec73313/calib/CsPad2x2::CalibV1/MecTargetChamber.0:Cspad2x2.1/'
    path  = '/reg/d/psdm/xpp/xpptut13/calib/CsPad2x2::CalibV1/XppGon.0:Cspad2x2.1/'
    run   = 123
    calib = CSPAD2x2CalibPars(path, run)
    coord = CSPAD2x2PixCoords(calib)
    coord.print_cspad2x2_geometry_pars()
    return coord, calib

#------------------------------

def test_0() :
    """Test of default constructor.
    """
    coord = CSPAD2x2PixCoords() 
    test_of_coord_arrs(coord)

#------------------------------

def test_1() :
    """Test of instantiation with external parameters.
    """
    coord = test_instantiation_1() 
    test_of_coord_arrs(coord)

#------------------------------

def test_2() :
    """Test of instantiation with calib=CSPAD2x2CalibPars(path, run).
    """
    coord, calib = test_instantiation_2() 
    test_of_coord_arrs(coord, calib)

#------------------------------

def test_3() :
    """Test of instantiation with external parameters.
    """
    coord = test_instantiation_1() 
    img2d = coord.get_cspad2x2_image(None)
    print 'img2d.shape =', img2d.shape
    
    gg.plotImageLarge(img2d, amp_range=(-1, 2), figsize=(12,11))
    gg.show()
 
#------------------------------

def test_4() :
    """Test of instantiation with calib=CSPAD2x2CalibPars(path, run).
    """
    coord, calib = test_instantiation_2() 
    test_of_image(coord, calib)

#------------------------------
 
if __name__ == "__main__" :
    if len(sys.argv)==1   : print 'Use command: python', sys.argv[0], '<test-number=0-5>'
    elif sys.argv[1]=='0' : test_0()
    elif sys.argv[1]=='1' : test_1()
    elif sys.argv[1]=='2' : test_2()
    elif sys.argv[1]=='3' : test_3()
    elif sys.argv[1]=='4' : test_4()
    elif sys.argv[1]=='5' : test_reshaping_arrs_for_cspad2x2()
    else : print 'Non-expected arguments: sys.argv=', sys.argv

    sys.exit ( 'End of test.' )

#------------------------------
