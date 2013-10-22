#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module CSPadArrayImage ...
#         Global methods for different style of presentation of cspad array as an image.
#
#------------------------------------------------------------------------

"""This module provides access to the calibration parameters

This software was developed for the SIT project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@see RelatedModule

@version $Id: 2013-08-19$

@author Mikhail S. Dubrovin
"""

#------------------------------
#  Module's version from CVS --
#------------------------------
__version__ = "$Revision: 4 $"
# $Source$

#------------------------------

import sys
import numpy as np

#------------------------------

def getCSPadArrayWithGap(arr, gap=3) :
    """In the CSPAD-size array (32*185,388) inserts the gap and returns the array with shape (4,8*185,388+gap)"""
    #print 'getCSPadArrayWithGap(...): Input array shape =', arr.shape
    nrows,ncols = arr.shape # (32*185,388) # <== expected input array shape
    if ncols != 388 or nrows<185 :
        print 'getCSPadArrayWithGap(...): WARNING! UNEXPECTED INPUT ARRAY SHAPE =', arr.shape
        return arr
    arr_gap = np.zeros( (nrows,gap), dtype=np.int16 )
    arr_halfs = np.hsplit(arr,2)
    arr_with_gap = np.hstack((arr_halfs[0], arr_gap, arr_halfs[1]))
    arr_with_gap.shape = (nrows,ncols+gap)
    return arr_with_gap


def getCSPadSegments2D(arr, gap=3, hspace=0) :
    """Returns the CSPAD array image of shape (8*185, 4*(388+gap)+3*hspace) with horizontal gaps and spaces"""
    arr_all = getCSPadArrayWithGap(arr, gap)
    arr_all.shape = (4,8*185,388+gap) # Reshape for quad index
    arr_hsp = np.zeros( (8*185, hspace), dtype=np.int16 )
    return np.hstack((arr_all[0,:],arr_hsp,arr_all[1,:],arr_hsp,arr_all[2,:],arr_hsp,arr_all[3,:]))


def getCSPadArrayAs2DImage(arr, gap=3, hspace=5, vspace=5) :
    """Returns the CSPAD array image of shape (8*185+7*vspace, 4*(388+gap)+3*hspace) with horizontal and vertical spaces"""
    arr_h = getCSPadSegments2D(arr, gap, hspace)
    hdim = (388+gap)*4 + hspace*3
    #print 'hdim =', hdim
    #print 'arr_h.shape =', arr_h.shape
    nrows,ncols = arr_h.shape
    arr_h.shape = (8, 185, ncols) # split for 2x1 index
    arr_vsp = np.zeros( (vspace, ncols), dtype=np.int16 )
    return np.vstack((arr_h[0,:], arr_vsp,\
                      arr_h[1,:], arr_vsp,\
                      arr_h[2,:], arr_vsp,\
                      arr_h[3,:], arr_vsp,\
                      arr_h[4,:], arr_vsp,\
                      arr_h[5,:], arr_vsp,\
                      arr_h[6,:], arr_vsp,\
                      arr_h[7,:]))

#------------------------------

import GlobalGraphics as gg # For test purpose in main only

#------------------------------

def get_raw_array_for_cspad_test() :
    """Returns raw cspad array for test purpose"""
    #arr = getRandomImage()
    arr = np.arange(32*185*388)
    arr.shape = (32*185, 388)
    return arr


def plot_img(img_arr) :
    """Plot image for test purpose"""
    print 'img_arr.shape=', img_arr.shape
    gg.plotImageLarge(img_arr) #, img_range=None, amp_range=None, ... 
    gg.savefig('cspad-arr-img.png')
    gg.move(500,10)
    gg.show()


def test1() :
    """Default pars for gaps"""
    img = getCSPadArrayAs2DImage(get_raw_array_for_cspad_test())
    plot_img(img)


def test2() :
    """Well-distinguished gaps"""
    img = getCSPadArrayAs2DImage(get_raw_array_for_cspad_test(),3,20,20)
    plot_img(img)


def test3() :
    """Test image with 0-gaps"""
    img = getCSPadArrayAs2DImage(get_raw_array_for_cspad_test(),0,0,0)
    plot_img(img)


def test4() :
    """Test helper getCSPadSegments2D"""
    img = getCSPadSegments2D(get_raw_array_for_cspad_test())
    plot_img(img)


def test5() :
    """Test helper getCSPadArrayWithGap"""
    img = getCSPadArrayWithGap(get_raw_array_for_cspad_test())
    plot_img(img)

#------------------------------

def main() :

    if len(sys.argv)==1   :
        print 'Use command > python %s <test-number [1-5]>' % sys.argv[0]
        sys.exit ('Add <test-number> in command line...')

    elif sys.argv[1]=='1' : test1()
    elif sys.argv[1]=='2' : test2()
    elif sys.argv[1]=='3' : test3()
    elif sys.argv[1]=='4' : test4()
    elif sys.argv[1]=='5' : test5()
    else :
        print 'Non-expected arguments: sys.argv=', sys.argv
        sys.exit ('Check input parameters')

#------------------------------

if __name__ == "__main__" :

    main()
    sys.exit ( 'End of test.' )

#------------------------------
