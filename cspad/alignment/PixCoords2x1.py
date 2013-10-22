#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module PixCoords2x1...
#
#------------------------------------------------------------------------

"""
This software was developed for the SIT project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@see RelatedModule

@version $Id: 2013-03-08$

@author Mikhail S. Dubrovin

Use matrix notations (like in data array)
DIFFERENT from the detector map... rows<->cols:
Assume that 2x1 has 195 rows and 388 columns
The (r,c)=(0,0) is in the top left corner of the matrix, has coordinates (xmin,ymax)

                    ^ Y          (Xmax,Ymax)
   (0,0)            |            (0,387)
      ------------------------------
      |             |              |
      |             |              |
      |             |              |
    --|-------------+--------------|----> X
      |             |              |
      |             |              |
      |             |              |
      ------------------------------
   (184,0)          |           (184,387)
   (Xmin,Ymin)

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
import GlobalGraphics as gg # For test purpose in main only
#------------------------------

def rotation(X, Y, C, S) :
    """For numpy arrays X and Y returns the numpy arrays of Xrot and Yrot
    """
    Xrot = X*C - Y*S 
    Yrot = Y*C + X*S 
    return Xrot, Yrot

#------------------------------

class PixCoords2x1() :
    """Self-sufficient class for generation of CSPad 2x1 sensor pixel coordinate array"""

    rows  = 185    # Number of rows in 2x1 at rotation 0
    cols  = 388    # Number of cols in 2x1 at rotation 0
    pixs  = 109.92 # Pixel size in um (micrometer)
    pixw  = 274.80 # Wide pixel size in um (micrometer)

    colsh = cols/2
    pixsh = pixs/2
    pixwh = pixw/2

#------------------------------

    def __init__(sp, use_wide_pix_center=True) :
        #print 'PixCoords2x1.__init__()'

        sp.use_wide_pix_center = use_wide_pix_center

        sp.x_map2x1_um_offset  = None
        sp.x_map2x1_pix_offset = None

        sp.make_maps_of_2x1_pix_coordinates()

#------------------------------

    def make_maps_of_2x1_pix_coordinates(sp) :
        """Makes [185,388] maps of x and y 2x1 pixel coordinates
        with origin in the center of 2x1
        """        
        x_rhs = np.arange(sp.colsh)*sp.pixs + sp.pixw - sp.pixsh
        if sp.use_wide_pix_center : x_rhs[0] = sp.pixwh # set x-coordinate of the wide pixel in its geometry center
        sp.x_arr_um = np.hstack([-x_rhs[::-1],x_rhs])

        sp.y_arr_um = -np.arange(sp.rows) * sp.pixs
        sp.y_arr_um -= sp.y_arr_um[-1]/2 # move origin to the center of array

        sp.x_arr_pix = sp.x_arr_um/sp.pixs
        sp.y_arr_pix = sp.y_arr_um/sp.pixs

        #sp.x_arr_pix = (sp.x_arr_um/sp.pixs + 0.25).astype(int)
        #sp.y_arr_pix = (sp.y_arr_um/sp.pixs + 0.25).astype(int)

        sp.x_map2x1_um,  sp.y_map2x1_um  = np.meshgrid(sp.x_arr_um,  sp.y_arr_um)
        sp.x_map2x1_pix, sp.y_map2x1_pix = np.meshgrid(sp.x_arr_pix, sp.y_arr_pix)
        
#------------------------------

    def print_maps_2x1_um(sp) :
        print 'x_map2x1_um = ',       sp.x_map2x1_um
        print 'x_map2x1_um.shape = ', sp.x_map2x1_um.shape
        print 'y_map2x1_um = ',       sp.y_map2x1_um
        print 'y_map2x1_um.shape = ', sp.y_map2x1_um.shape

#------------------------------

    def print_maps_2x1_pix(sp) :
        print 'x_map2x1_pix = ',       sp.x_map2x1_pix
        print 'x_map2x1_pix.shape = ', sp.x_map2x1_pix.shape
        print 'y_map2x1_pix = ',       sp.y_map2x1_pix
        print 'y_map2x1_pix.shape = ', sp.y_map2x1_pix.shape

#------------------------------

    def print_xy_arr_um(sp) :
        print 'x_arr_um:\n',       sp.x_arr_um
        print 'x_arr_um.shape = ', sp.x_arr_um.shape
        print 'y_arr_um:\n',       sp.y_arr_um
        print 'y_arr_um.shape = ', sp.y_arr_um.shape

#------------------------------

    def print_xy_arr_pix(sp) :
        print 'x_arr_pix:\n',       sp.x_arr_pix
        print 'x_arr_pix.shape = ', sp.x_arr_pix.shape
        print 'y_arr_pix:\n',       sp.y_arr_pix
        print 'y_arr_pix.shape = ', sp.y_arr_pix.shape

#------------------------------

    def print_xy_min_max_um(sp) :
        xmin, ymin = sp.get_xy_min_um()
        xmax, ymax = sp.get_xy_max_um()
        print 'In [um] xmin:%9.2f, xmax:%9.2f, ymin:%9.2f, ymax:%9.2f' % (xmin, xmax, ymin, ymax)

#------------------------------

    def print_xy_min_max_pix(sp) :
        xmin, ymin = sp.get_xy_min_pix()
        xmax, ymax = sp.get_xy_max_pix()
        print 'In [pix] xmin:%5.0f, xmax:%5.0f, ymin:%5.0f, ymax:%5.0f' % (xmin, xmax, ymin, ymax)

#------------------------------

    def get_xy_min_um(sp) : 
        return sp.x_arr_um[0], sp.y_arr_um[-1]

    def get_xy_max_um(sp) : 
        return sp.x_arr_um[-1], sp.y_arr_um[0]

    def get_xy_min_pix(sp) : 
        return sp.x_arr_pix[0], sp.y_arr_pix[-1]

    def get_xy_max_pix(sp) : 
        return sp.x_arr_pix[-1], sp.y_arr_pix[0]

    def get_cspad2x1_xy_maps_um(sp) : 
        return sp.x_map2x1_um, sp.y_map2x1_um

    def get_cspad2x1_xy_maps_pix(sp) : 
        return sp.x_map2x1_pix, sp.y_map2x1_pix

    def get_cspad2x1_xy_maps_um_with_offset(sp) : 
        if  sp.x_map2x1_um_offset == None :
            x_min_um, y_min_um = sp.get_xy_min_um()
            sp.x_map2x1_um_offset = sp.x_map2x1_um - x_min_um
            sp.y_map2x1_um_offset = sp.y_map2x1_um - y_min_um
        return sp.x_map2x1_um_offset, sp.y_map2x1_um_offset

    def get_cspad2x1_xy_maps_pix_with_offset(sp) : 
        if  sp.x_map2x1_pix_offset == None :
            x_min_pix, y_min_pix = sp.get_xy_min_pix()
            sp.x_map2x1_pix_offset = sp.x_map2x1_pix - x_min_pix
            sp.y_map2x1_pix_offset = sp.y_map2x1_pix - y_min_pix
        return sp.x_map2x1_pix_offset, sp.y_map2x1_pix_offset

#------------------------------
# cspad2x1 = PixCoords2x1(use_wide_pix_center=False)
#------------------------------
#------------------------------
#------------------------------
#----------- TEST -------------
#------------------------------
#------------------------------
#------------------------------

def test_2x1_xy_maps() :

    w = PixCoords2x1()
    w.print_maps_2x1_um()

    titles = ['X map','Y map']
    #for i,arr2d in enumerate([w.x_map2x1,w.y_map2x1]) :
    for i,arr2d in enumerate( w.get_cspad2x1_xy_maps_pix() ) :
        amp_range = (arr2d.min(), arr2d.max())
        gg.plotImageLarge(arr2d, amp_range=amp_range, figsize=(10,5), title=titles[i])
        gg.move(200*i,100*i)

    gg.show()

#------------------------------

def test_2x1_img() :

    t0_sec = time()
    w = PixCoords2x1(use_wide_pix_center=False)
    #w = PixCoords2x1(use_wide_pix_center=True)
    print 'Consumed time for coordinate arrays (sec) =', time()-t0_sec

    X,Y = w.get_cspad2x1_xy_maps_pix()
    w.print_xy_arr_um()
    w.print_xy_arr_pix()
    w.print_xy_min_max_um()
    w.print_xy_min_max_pix()


    #print 'X(pix) :\n', X
    print 'X.shape =', X.shape

    xmin, ymin = w.get_xy_min_pix()
    xmax, ymax = w.get_xy_max_pix()
    xmin-=0.5; xmax+=0.5; ymin-=0.5; ymax+=0.5;

    xsize = xmax - xmin 
    ysize = ymax - ymin 
    print 'xsize =', xsize # 391.0 
    print 'ysize =', ysize # 185.0

    H, Xedges, Yedges = np.histogram2d(X.flatten(), Y.flatten(), bins=[xsize,ysize], range=[[xmin, xmax], [ymin, ymax]], normed=False, weights=X.flatten()+Y.flatten()) 

    print 'Xedges:', Xedges
    print 'Yedges:', Yedges
    print 'H.shape:', H.shape

    gg.plotImageLarge(H, amp_range=(-250, 250), figsize=(8,10)) # range=(-1, 2), 
    gg.show()

#------------------------------

def test_2x1_img_easy() :
    pc2x1 = PixCoords2x1(use_wide_pix_center=False)
    #X,Y = pc2x1.get_cspad2x1_xy_maps_pix()
    X,Y = pc2x1.get_cspad2x1_xy_maps_pix_with_offset()
    iX, iY = (X+0.25).astype(int), (Y+0.25).astype(int)
    img = gg.getImageFromIndexArrays(iX,iY,iX+iY)
    gg.plotImageLarge(img, amp_range=(0, 500), figsize=(8,10))
    gg.show()
 
#------------------------------
 
if __name__ == "__main__" :

    if len(sys.argv)==1   : print 'For other test(s) use command: python', sys.argv[0], '<test-number=1-3>'
    elif sys.argv[1]=='1' : test_2x1_xy_maps()
    elif sys.argv[1]=='2' : test_2x1_img()
    elif sys.argv[1]=='3' : test_2x1_img_easy()
    else : print 'Non-expected arguments: sys.argv=', sys.argv

    sys.exit( 'End of test.' )

#------------------------------
