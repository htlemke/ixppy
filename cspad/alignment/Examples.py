#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module Examples...
#
#------------------------------------------------------------------------

"""This module provides examples of how to get and use the CSPad image

This software was developed for the SIT project.  If you use all or
part of it, please give an appropriate acknowledgment.

@see RelatedModule

@version $Id: 2011-11-18$

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
import numpy              as np

import PyCSPadImage.CalibPars          as calp
import PyCSPadImage.CSPadConfigPars    as ccp
import PyCSPadImage.CSPadImageProducer as cip
import PyCSPadImage.CSPADPixCoords     as pixcoor
import PyCSPadImage.PixCoords2x1       as pixcoor2x1


import PyCSPadImage.GlobalGraphics     as gg # For test purpose in main only
import PyCSPadImage.GlobalMethods      as gm # For test purpose in main only
import PyCSPadImage.HDF5Methods        as hm # For test purpose in main only

#------------------------------
#------------------------------
#------------------------------
#------------------------------

def test_plot_cspad_image(fname, dsname, path_calib, run, event=0, nevents=1, amps=None, do_peds=False) :
    """Test of instantiation with external parameters.
    """
    
    calib  = calp.CalibPars(path_calib, run)
    coord  = pixcoor.CSPADPixCoords(calib)
    coord.print_cspad_geometry_pars()

    ds1ev  = None
    if nevents == 1 : ds1ev = hm.getOneCSPadEventForTest( fname, dsname, event )
    else            : ds1ev = hm.getAverageCSPadEvent   ( fname, dsname, event1=event, nevents=nevents )

    if do_peds :
        peds = calib.getCalibPars('pedestals')
        peds.shape = (32, 185, 388)
        ds1ev -= peds

        #ped_fname = '/reg/neh/home1/dubrovin/LCLS/CSPadPedestals/cspad-pedestals-xpp66213-r0149.dat' # shape = (5920, 388)
        #peds = gm.getCSPadArrayFromFile(ped_fname)
        #print 'peds.shape:', peds.shape # (32, 185, 388)

    config = ccp.CSPadConfigPars()
    config.setCSPadConfiguration( fname, dsname, event ) # This will set current CSPad configuration
    config.printCSPadConfigPars()
 
    img2d = coord.get_cspad_image(ds1ev, config)
    print 'img2d.shape =', img2d.shape
    
    gg.plotImageLarge(img2d, amp_range=amps, figsize=(12,11))
    #gg.plotImageLarge(img2d, amp_range=None, figsize=(12,11))
    gg.savefig('cspad-img.png')
    gg.show()

#------------------------------

def test_cspad_image(test_num=1) :
    """Test of instantiation with external parameters.
    """

    event       = 10
    nevents     = 1
    dsname      = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/CxiDs1.0:Cspad.0/data'
    do_peds     = False

    if test_num == 1 : # Wather ring and shadows 
        run        = 150
        fname      = '/reg/d/psdm/xpp/xpptut13/hdf5/xpptut13-r0150.h5' # xpp66213-r0150.h5
        path_calib = '/reg/d/psdm/xpp/xpptut13/calib/CsPad::CalibV1/XppGon.0:Cspad.0/'
        amps       = (-10, 200)
        # XPP specific:
        do_peds    = True
        dsname     = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/XppGon.0:Cspad.0/data'

    elif test_num == 2 : # Wires
        run        = 9
        fname      = '/reg/d/psdm/cxi/cxitut13/hdf5/cxitut13-r0009.h5' # cxi35711-r0009.h5
        path_calib = '/reg/d/psdm/cxi/cxi35711/calib/CsPad::CalibV1/CxiDs1.0:Cspad.0/'
        amps       = (0, 2000)

    elif test_num == 3 : # Wires and missing 2x1 
        run        = 628
        fname      = '/reg/d/psdm/cxi/cxitut13/hdf5/cxitut13-r0628.h5' # cxi80410-r0628.h5
        path_calib = '/reg/d/psdm/cxi/cxi80410/calib/CsPad::CalibV1/CxiDs1.0:Cspad.0/'
        amps       = (0, 3000)
    
    elif test_num == 4 : # Equidistant rings
        run        = 1150
        fname      = '/reg/d/psdm/cxi/cxitut13/hdf5/cxitut13-r1150.h5' # cxi80410-r1150.h5
        path_calib = '/reg/d/psdm/cxi/cxi80410/calib/CsPad::CalibV1/CxiDs1.0:Cspad.0/'
        amps       = (0, 200)

    elif test_num == 5 : # Rings ? (T.J.) 
        run        = 135
        fname      = '/reg/d/psdm/cxi/cxitut13/hdf5/cxitut13-r0135.h5' # cxi64813-r0135.h5
        path_calib = '/reg/d/psdm/cxi/cxi64813/calib/CsPad::CalibV1/CxiDs1.0:Cspad.0/'
        amps       = (-10, 100)

    elif test_num == 6 : # Test of T.J. alignment
        print 'HERE!'
        run        = 13
        fname      = '/reg/d/psdm/CXI/cxia4113/hdf5/cxia4113-r0013.h5'
        path_calib = '/reg/neh/home1/dubrovin/LCLS/CSPadAlignment-v01/calib-test-cxia4113-r0013-Ds1/CsPad::CalibV1/CxiDs1.0:Cspad.0/'
        do_peds    = True
        event      = 600
        nevents    = 100
        amps       = (0, 500)

    else: 
        print 'Non-defined test number:', test_num
        sys.exit ( 'Exit, use proper input parameter.' )        

    test_plot_cspad_image(fname, dsname, path_calib, run, event, nevents, amps, do_peds)

#------------------------------

if __name__ == "__main__" :
    if len(sys.argv)==1   :
        print 'Use command: python', sys.argv[0], '<test-number=1-5, 10-13, 20-22>'
        sys.exit ( 'Exit, use proper input parameter.' )        

    test_number = int(sys.argv[1])

    if   test_number < 10 : test_cspad_image(test_number)

    elif sys.argv[1]=='10': pixcoor.test_cspadpixcoords_0()
    elif sys.argv[1]=='11': pixcoor.test_cspadpixcoords_1()
    elif sys.argv[1]=='12': pixcoor.test_cspadpixcoords_2()
    elif sys.argv[1]=='13': pixcoor.test_cspadpixcoords_3()

    elif sys.argv[1]=='20': pixcoor2x1.test_2x1_xy_maps()
    elif sys.argv[1]=='21': pixcoor2x1.test_2x1_img()
    elif sys.argv[1]=='22': pixcoor2x1.test_2x1_img_easy()()

    else : print 'Non-expected arguments: sys.argv=', sys.argv

    sys.exit ( 'End of test.' )

#----------------------------------------------
