#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module CSPadPixCoordsWODB...
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

#import matplotlib.pyplot as plt
import GlobalGraphics     as gg # For test purpose in main only
#------------------------------

def rotation(X, Y, C, S) :
    """For numpy arrays X and Y returns the numpy arrays of Xrot and Yrot
    """
    Xrot = X*C + Y*S 
    Yrot = Y*C - X*S 
    return Xrot, Yrot

#------------------------------

class CSPadPixCoordsWODB () :
    """Self-sufficient class for generation of CSPad pixel coordinate array without data base (WODB)"""

    quads       =   4    # Total number of quads in cspad
    sects       =   8    # Total number of sections in quad
    rows        = 185    # Number of rows in 2x1 at rotation 0
    cols        = 388    # Number of cols in 2x1 at rotation 0
    pixs        = 109.92 # Pixel size in um (micrometer)
    pixw        = 274.80 # Wide pixel size in um (micrometer)

    colsh = cols/2
    pixsh = pixs/2
    pixwh = pixw/2

#------------------------------

    def __init__ (sp, xc_um=None, yc_um=None, orient_deg=None, tilt_deg=None) :
        print 'CSPadPixCoordsWODB __init__'

        if xc_um      == None : return
        if yc_um      == None : return
        if orient_deg == None : return
        if tilt_deg   == None : return

        sp.make_cspad_pix_coordinate_arrays (xc_um, yc_um, orient_deg, tilt_deg)

#------------------------------

    def make_maps_of_2x1_pix_coordinates (sp) :
        """Makes [185,388] maps of x and y 2x1 pixel coordinates"""        
        x_rhs = np.arange(sp.colsh)*sp.pixs + sp.pixw - sp.pixsh
        x_rhs[0] = sp.pixwh # set x-coordinate of the wide pixel 
        x_arr = np.hstack([-x_rhs[::-1],x_rhs])

        y_arr = np.arange(sp.rows) * sp.pixs
        y_arr -= y_arr[-1]/2 # move origin to the center of array

        sp.x_map2x1, sp.y_map2x1 = np.meshgrid(x_arr, y_arr)

#------------------------------

    def print_maps_2x1(sp) :
        print 'x_map2x1 = ',       sp.x_map2x1
        print 'x_map2x1.shape = ', sp.x_map2x1.shape
        print 'y_map2x1 = ',       sp.y_map2x1
        print 'y_map2x1.shape = ', sp.y_map2x1.shape

#------------------------------

    def make_cspad_pix_coordinate_arrays (sp, xc_um, yc_um, orient_deg, tilt_deg) : # All lists of [4,8]
        """Makes [4,8,185,388] cspad pixel x and y coordinate arrays"""        
        sp.make_maps_of_2x1_pix_coordinates()

        sp.x_pix_um = np.zeros((sp.quads,sp.sects,sp.rows,sp.cols), dtype=np.float32)
        sp.y_pix_um = np.zeros((sp.quads,sp.sects,sp.rows,sp.cols), dtype=np.float32)

        angle_deg = orient_deg + tilt_deg
 
        for quad in range(sp.quads) :
            for sect in range(sp.sects) :

                angle_rad = math.radians(angle_deg[quad][sect])                
                S,C = math.sin(angle_rad), math.cos(angle_rad)
                Xrot, Yrot = rotation(sp.x_map2x1, sp.y_map2x1, C, S)

                sp.x_pix_um[quad][sect][:] =  Xrot + xc_um[quad][sect]
                sp.y_pix_um[quad][sect][:] =  Yrot + yc_um[quad][sect]

        sp.x_pix_um -= sp.x_pix_um.min() + 5 # add offset in um to get rid of "rounding" strips...
        sp.y_pix_um -= sp.y_pix_um.min() + 5

#------------------------------

    def get_cspad_pix_coordinate_arrays_um (sp) : 
        return sp.x_pix_um, sp.y_pix_um


    def get_cspad_pix_coordinate_arrays_pix (sp) : 
        return sp.x_pix_um/sp.pixs, sp.y_pix_um/sp.pixs

#------------------------------

    def print_cspad_coordinate_arrays(sp) :
        print 'sp.x_pix_um:\n',      sp.x_pix_um
        print 'sp.x_pix_um.shape =', sp.x_pix_um.shape
        print 'sp.y_pix_um\n',       sp.y_pix_um
        print 'sp.y_pix_um.shape =', sp.y_pix_um.shape

#------------------------------
#------------------------------
#------------------------------
#----------- TEST -------------
#------------------------------
#------------------------------
#------------------------------

def main_test_2x1() :

    w = CSPadPixCoordsWODB()
    w.make_maps_of_2x1_pix_coordinates()
    w.print_maps_2x1()

    for i,arr2d in enumerate([w.x_map2x1,w.y_map2x1]) :
        range = (arr2d.min(), arr2d.max())
        gg.plotImage(arr2d, range, figsize=(10,5))
        gg.move(200*i,100*i)

    gg.show()

#------------------------------

def main_test_cspad() :

    xc_um = np.array(
            [[ 473.38,  685.26,  155.01,  154.08,  266.81,   53.95,  583.04,  582.15],  
             [ 989.30,  987.12, 1096.93,  884.11, 1413.16, 1414.94, 1500.83, 1288.02],  
             [1142.59,  930.23, 1459.44, 1460.67, 1347.57, 1559.93, 1032.27, 1033.44],  
             [ 626.78,  627.42,  516.03,  729.15,  198.28,  198.01,  115.31,  327.66]]) * 109.92

    yc_um = np.array(
            [[1028.07, 1026.28, 1139.46,  926.91, 1456.78, 1457.35, 1539.71, 1327.89],  
             [1180.51,  967.36, 1497.74, 1498.54, 1385.08, 1598.19, 1069.65, 1069.93],  
             [ 664.89,  666.83,  553.60,  765.91,  237.53,  236.06,  152.17,  365.47],  
             [ 510.38,  722.95,  193.33,  193.41,  308.04,   95.25,  625.28,  624.14]]) * 109.92

    orient_deg = np.array(
                    [[  90.,   90.,    0.,    0.,  270.,  270.,    0.,    0.],
                     [   0.,    0.,  270.,  270.,  180.,  180.,  270.,  270.],
                     [  90.,   90.,    0.,    0.,  270.,  270.,    0.,    0.],
                     [   0.,    0.,  270.,  270.,  180.,  180.,  270.,  270.]])
 
    tilt_deg = np.array(
                    [[0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],  
                     [0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],  
                     [0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],  
                     [0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]])

    print 'xc_um:\n',      xc_um
    print 'yc_um:\n',      yc_um
    print 'orient_deg:\n', orient_deg
    print 'tilt_deg:\n',   tilt_deg

    w = CSPadPixCoordsWODB(xc_um, yc_um, orient_deg, tilt_deg)
    #w.make_cspad_pix_coordinate_arrays (xc_um, yc_um, orient_deg, tilt_deg)
    w.print_cspad_coordinate_arrays()
    X,Y = w.get_cspad_pix_coordinate_arrays_pix ()

    print 'X(pix) :\n', X
    print 'X.shape =\n', X.shape

    xsize = X.max() + 1
    ysize = Y.max() + 1
    H,Xedges,Yedges = np.histogram2d(X.flatten(), Y.flatten(), bins=[xsize,ysize], range=[[0,xsize],[0,ysize]], normed=False, weights=None) 

    range = (-5, 5)
    gg.plotImageLarge(H, range=(-1, 2), figsize=(12,11))
    gg.show()

#------------------------------
 
if __name__ == "__main__" :

    #main_test_2x1()

    main_test_cspad()

    sys.exit ( 'End of test.' )

#------------------------------
