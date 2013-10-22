#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module CSPADPixCoords...
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
import GlobalGraphics  as gg # For test purpose in main only
#------------------------------

class CSPADPixCoords (PixCoords2x1) :
    """Class for generation of CSPad pixel coordinate array with and without data base

       Interface
       =========
       1.0 Instantiation with default parameters taken from optical measurement for XPP on 2013-01-29:
           coord = CSPADPixCoords()

       1.1 Instantiation with external geometry parameters:
 
           All parameters optional. Default values will be used if parameters are not specified.
           xc    - np.array(...), shape=(3, 4, 8)) [um]
           yc    - np.array(...), shape=(3, 4, 8)) [um]
           tilt  - np.array(...), shape=(4, 8)) [deg]
           coord = CSPADPixCoords(xc_um=xc, yc_um=yc, tilt_deg=tilt)


       1.2 Instantiation with regular calibration parameters:

           path  = '/reg/neh/home1/dubrovin/LCLS/CSPadAlignment-v01/calib-xpp-2013-01-29'
           run   = 123
           calib = CalibPars(path, run)
           coord = CSPADPixCoords(calib)

       1.2 Access methods:
           Get arrays of pixel coordinates in mu with shape: [2,185,388]
            X, Y = coord.get_cspad_pix_coordinate_arrays_um  (config=None)
           or in integer pixels:
           iX,iY = coord.get_cspad_pix_coordinate_arrays_pix (config=None)

       1.3 Get image  
           img = coord.get_cspad_image(data,config)       
    """

    quads = 4 # Total number of quads in cspad
    sects = 8 # Total number of sections in quad

    # Default from optical measurement for XPP on 2013-01-29
    xc_um_def  = np.array(
          [[ 477.78,    690.20,    159.77,    160.06,    277.17,     64.77,    591.30,    591.01],
           [ 990.78,    989.30,   1105.38,    891.19,   1421.65,   1423.66,   1502.28,   1289.93],
           [1143.85,    932.00,   1461.86,   1463.74,   1349.75,   1562.62,   1032.39,   1033.60],
           [ 633.06,    632.80,    518.88,    731.75,    200.62,    198.75,    118.50,    331.23]] ) * PixCoords2x1.pixs # 109.92

    yc_um_def = np.array(
          [[1018.54,   1019.42,   1134.27,    921.94,   1451.06,   1451.01,   1532.55,   1319.23],
           [1173.24,    960.71,   1490.18,   1491.45,   1374.97,   1587.78,   1058.56,   1061.14],
           [ 658.23,    658.54,    542.73,    755.26,    225.91,    224.22,    146.39,    358.27],
           [ 507.44,    720.59,    189.73,    190.28,    306.25,     93.65,    620.68,    619.85]] ) * PixCoords2x1.pixs # 109.92

    zc_um_def = np.zeros((4,8), dtype=np.float32)

    # Orientation of sect:   0   1     2     3     4     5     6     7
    orient_opt = np.array ( [0., 0., 270., 270., 180., 180., 270., 270.] )

     # Orientation of quad:         0              1           2              3 
    orient_def = np.array ( [orient_opt+90, orient_opt, orient_opt-90, orient_opt-180] )

    # Tilt:
    tilt_2013_01_29 = np.array (
          [[  0.27766,   0.37506,   0.11976,   0.17369,  -0.04934,   0.01119,   0.13752,  -0.00921],   
           [ -0.45066,  -0.18880,  -0.20400,  -0.33507,  -0.62242,  -0.40196,  -0.56593,  -0.59475],   
           [ -0.03290,   0.00658,  -0.33954,  -0.27106,  -0.71923,  -0.31647,   0.02829,   0.10723],   
           [ -0.11054,   0.10658,   0.25005,   0.16121,  -0.58560,  -0.43369,  -0.26916,  -0.18225]] )

    #tilt_def = tilt_2013_01_29
    tilt_def = np.zeros((quads, sects), dtype=np.float32)

#------------------------------

    def __init__ (sp, calib=None, xc_um=None, yc_um=None, tilt_deg=None, use_wide_pix_center=False) :
        #print 'CSPAD2x2PixCoords.__init__(...)'

        PixCoords2x1.__init__ (sp, use_wide_pix_center)

        if calib == None :

            if xc_um == None : sp.xc = sp.xc_um_def
            else             : sp.xc = xc_um

            if yc_um == None : sp.yc = sp.yc_um_def
            else             : sp.yc = yc_um

            sp.zc   = sp.zc_um_def 

            if tilt_deg == None : sp.tilt = sp.tilt_def
            else                : sp.tilt = tilt_deg 

        else :
            [sp.xc,sp.yc,sp.zc] = calib.getCalibPars('center_global') * PixCoords2x1.pixs # 109.92
            sp.tilt             = calib.getCalibPars('tilt')
            #print 'USE xc, yc, zc =', sp.xc, sp.yc, sp.zc

        sp.orient = sp.orient_def
        sp.calib  = calib

        sp.make_cspad_pix_coordinate_arrays (sp.xc, sp.yc, sp.orient, sp.tilt)

#------------------------------

    def print_cspad_geometry_pars (sp) :
        msg = 'print_cspad_geometry_pars():' \
            + '\nxc [pix]:\n' + str( sp.xc/PixCoords2x1.pixs ) \
            + '\nyc [pix]:\n' + str( sp.yc/PixCoords2x1.pixs ) \
            + '\norient:\n'   + str( sp.orient ) \
            + '\ntilt:\n'     + str( sp.tilt )
            #+ '\nxc:'       + str( sp.xc ) \
            #+ '\nyc:'       + str( sp.yc ) \
        print msg

#------------------------------

    def make_cspad_pix_coordinate_arrays (sp, xc_um, yc_um, orient_deg, tilt_deg=None) : # All lists of [4,8]
        """Makes [4,8,185,388] cspad pixel x and y coordinate arrays"""        
        sp.make_maps_of_2x1_pix_coordinates()

        sp.x_pix_um = np.zeros((sp.quads,sp.sects,sp.rows,sp.cols), dtype=np.float32)
        sp.y_pix_um = np.zeros((sp.quads,sp.sects,sp.rows,sp.cols), dtype=np.float32)

        angle_deg = np.array(orient_deg)
        if tilt_deg != None : angle_deg += tilt_deg
 
        for quad in range(sp.quads) :
            for sect in range(sp.sects) :

                angle_rad = math.radians(angle_deg[quad][sect])                
                S,C = math.sin(angle_rad), math.cos(angle_rad)
                Xrot, Yrot = rotation(sp.x_map2x1_um, sp.y_map2x1_um, C, S) # defined in PixCoords2x1

                sp.x_pix_um[quad][sect][:] =  Xrot + xc_um[quad][sect]
                sp.y_pix_um[quad][sect][:] =  Yrot + yc_um[quad][sect]

        sp.x_pix_um -= sp.x_pix_um.min()
        sp.y_pix_um -= sp.y_pix_um.min()

        sp.x_pix_pix = (sp.x_pix_um/sp.pixs+0.25).astype(int) 
        sp.y_pix_pix = (sp.y_pix_um/sp.pixs+0.25).astype(int)

#------------------------------

    def get_cspad_pix_arrays_shaped_by_config (sp, x_arr, y_arr, config=None) :
        """Array shaping as data in case if the config parameter is defined: config = CSPadConfigPars()..."""
        if config == None :
            return x_arr, y_arr
        else :
            sp.x_arr_as_data = config.getCSPadPixArrayShapedAsData(x_arr)
            sp.y_arr_as_data = config.getCSPadPixArrayShapedAsData(y_arr)       
            return sp.x_arr_as_data, sp.y_arr_as_data

#------------------------------

    def get_cspad_pix_coordinate_arrays_um (sp, config=None) : 
        return sp.get_cspad_pix_arrays_shaped_by_config (sp.x_pix_um, sp.y_pix_um, config)
    

    def get_cspad_pix_coordinate_arrays_pix (sp, config=None) : 
        return sp.get_cspad_pix_arrays_shaped_by_config (sp.x_pix_pix, sp.y_pix_pix, config)

#------------------------------

    def print_cspad_coordinate_arrays(sp) :
        print 'sp.x_pix_um:\n',      sp.x_pix_um
        print 'sp.x_pix_um.shape =', sp.x_pix_um.shape
        print 'sp.y_pix_um\n',       sp.y_pix_um
        print 'sp.y_pix_um.shape =', sp.y_pix_um.shape

#------------------------------

    def get_cspad_image(sp, data_arr=None, config=None) : # preferable data_arr.shape=(32,185,388) like in data
        """ Test of coordinate arrays, plot image map.
        """
        iX,iY = sp.get_cspad_pix_coordinate_arrays_pix(config)
        data = data_arr
        #if data_arr != None and data.shape != (4,8,185,388) :
        #    data.shape == (4,8,185,388)
        if data_arr != None : data = data_arr.flatten()
        return gg.getImageFromIndexArrays(iX.flatten(),iY.flatten(),data) # All arrays should have the same shape

#------------------------------
#------------------------------
#------------------------------
#----------- TEST -------------
#------------------------------
#------------------------------
#------------------------------

def get_test_cspad_pix_arr(config=None) :
    secs, rows, cols = shape = (32, 185, 388)
    arr_cspad_pix = np.zeros(shape, dtype=np.float32)
    arr_2x1_pix = np.ones((rows, cols), dtype=np.float32)
    for s in range(secs) :
        factor = s+4
        if s%2 : factor += 4
        arr_cspad_pix[s,:] += factor*arr_2x1_pix[:]

    arr_cspad_pix.shape = (4, 8, rows, cols)

    arr_out = arr_cspad_pix
    if config != None : arr_out = config.getCSPadPixArrayShapedAsData(arr_cspad_pix)
    return np.array(arr_out)

#------------------------------

def test_of_coord_arrs(coord, config=None) :
    """ Test of coordinate arrays, plot image map.
    """

    iX,iY = coord.get_cspad_pix_coordinate_arrays_pix (config)

    print 'iX.shape =', iX.shape
    print 'iY.shape =', iY.shape

    weights = get_test_cspad_pix_arr(config)

    t0_sec = time()
    #img2d = gg.getImageAs2DHist(iX,iY,W=None)
    img2d = gg.getImageFromIndexArrays(iX,iY,W=weights)
    print 'Consumed time to create image (sec) =', time()-t0_sec

    #gg.plotImageLarge(img2d, amp_range=(-1, 32), figsize=(12,11)) #amp_range=(0, 2000)
    gg.plotImageLarge(img2d, amp_range=None, figsize=(12,11)) #amp_range=(0, 2000)
    gg.show()

#------------------------------

def test_cspadpixcoords_instantiation_1() :
    """ Instantiation with external geometry parameters.
    """
    # Optical measurement for XPP from 2013-01-29
    xc = np.array(
          [[ 477.78,    690.20,    159.77,    160.06,    277.17,     64.77,    591.30,    591.01],
           [ 990.78,    989.30,   1105.38,    891.19,   1421.65,   1423.66,   1502.28,   1289.93],
           [1143.85,    932.00,   1461.86,   1463.74,   1349.75,   1562.62,   1032.39,   1033.60],
           [ 633.06,    632.80,    518.88,    731.75,    200.62,    198.75,    118.50,    331.23]] ) * PixCoords2x1.pixs # 109.92

    yc = np.array(
          [[1018.54,   1019.42,   1134.27,    921.94,   1451.06,   1451.01,   1532.55,   1319.23],
           [1173.24,    960.71,   1490.18,   1491.45,   1374.97,   1587.78,   1058.56,   1061.14],
           [ 658.23,    658.54,    542.73,    755.26,    225.91,    224.22,    146.39,    358.27],
           [ 507.44,    720.59,    189.73,    190.28,    306.25,     93.65,    620.68,    619.85]] ) * PixCoords2x1.pixs # 109.92

    tilt = np.array([[  0.27766,   0.37506,   0.11976,   0.17369,  -0.04934,   0.01119,   0.13752,  -0.00921],   
                     [ -0.45066,  -0.18880,  -0.20400,  -0.33507,  -0.62242,  -0.40196,  -0.56593,  -0.59475],   
                     [ -0.03290,   0.00658,  -0.33954,  -0.27106,  -0.71923,  -0.31647,   0.02829,   0.10723],   
                     [ -0.11054,   0.10658,   0.25005,   0.16121,  -0.58560,  -0.43369,  -0.26916,  -0.18225]] )

    #tilt = np.zeros((4, 8), dtype=np.float32)

    t0_sec = time()
    coord = CSPADPixCoords(xc_um=xc, yc_um=yc, tilt_deg=tilt, use_wide_pix_center=False)
    coord.print_cspad_geometry_pars()
    print 'Consumed time for CSPADPixCoords instatiation (sec) =', time()-t0_sec
    return coord

#------------------------------

import PyCSPadImage.CSPadConfigPars as ccp # For test purpose only
from   PyCSPadImage.CalibPars import *

def test_cspadpixcoords_instantiation_2() :
    """ Instantiation with regular calibration parameters.
    """
    path = '/reg/neh/home1/dubrovin/LCLS/CSPadAlignment-v01/calib-xpp-2013-01-29'
    run   = 123
    calib = CalibPars(path, run)
    print 'center_global:\n', calib.getCalibPars ('center_global') 
    coord = CSPADPixCoords(calib)
    coord.print_cspad_geometry_pars()
    return coord

#------------------------------

def test_cspadpixcoords_0() :
    """Test of default constructor.
    """
    coord = CSPADPixCoords() 
    coord.print_cspad_geometry_pars()
    test_of_coord_arrs(coord)

#------------------------------

def test_cspadpixcoords_1() :
    """Test of instantiation with external parameters.
    """
    coord = test_cspadpixcoords_instantiation_1() 
    test_of_coord_arrs(coord)

#------------------------------

def test_cspadpixcoords_2() :
    """Test of instantiation with calib=CSPADCalibPars(path, run).
    """
    coord = test_cspadpixcoords_instantiation_2() 
    coord.print_cspad_geometry_pars()
    test_of_coord_arrs(coord)

#------------------------------

def test_cspadpixcoords_3() :
    """Test of instantiation with external parameters.
    """
    coord = test_cspadpixcoords_instantiation_1() 
    img2d = coord.get_cspad_image(None)
    print 'img2d.shape =', img2d.shape
    
    gg.plotImageLarge(img2d, amp_range=(-1, 2), figsize=(12,11))
    gg.show()

#------------------------------

def test_cspadpixcoords_4() :
    """Test of default constructor for coords and non-default for config.
    """
    coord = CSPADPixCoords() 
    coord.print_cspad_geometry_pars()

    config = ccp.CSPadConfigPars() # instatiate object
    indPairs = np.array( 
        [ [ 0,   1,  -1,   3,   4,   5,   6,   7],
          [ 8,   9,  10,  11,  12,  13,  14,  15],
          [16,  17,  18,  19,  20,  21,  22,  23],
          [24,  -1,  25,  26,  27,  -1,  29,  30] ] )
    quadNums = [2, 3, 0, 1]
    config.setCSPadConfigArrays( indPairsInQuads=indPairs, quadNumsInEvent=quadNums )
    config.printCSPadConfigPars()

    test_of_coord_arrs(coord, config)

#------------------------------

if __name__ == "__main__" :
    if len(sys.argv)==1   : print 'Use command: python', sys.argv[0], '<test-number=0-4>'
    elif sys.argv[1]=='0' : test_cspadpixcoords_0() # Instatiation default
    elif sys.argv[1]=='1' : test_cspadpixcoords_1() # Instatiation using external geometry parameters
    elif sys.argv[1]=='2' : test_cspadpixcoords_2() # Instatiation using calib = CalibPars(path, run)
    elif sys.argv[1]=='3' : test_cspadpixcoords_3() # Test of coord.get_cspad_image()
    elif sys.argv[1]=='4' : test_cspadpixcoords_4() # Test of default constructor for coords and non-default for config
    else : print 'Non-expected arguments: sys.argv=', sys.argv

    sys.exit ( 'End of test.' )

#------------------------------
