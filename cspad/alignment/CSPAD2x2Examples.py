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

@version $Id: 2013-05-10$

@author Mikhail S. Dubrovin
"""

#--------------------------------
#  Module's version from CVS --
#--------------------------------
__version__ = "$Revision: 4 $"
# $Source$
#--------------------------------
import sys
import numpy as np

import PyCSPadImage.CSPAD2x2PixCoords        as pixcoor
import PyCSPadImage.CSPAD2x2CalibPars        as calpars
import PyCSPadImage.CSPAD2x2CalibParsDefault as cpd
import PyCSPadImage.PixCoords2x1             as pixcoor2x1

import PyCSPadImage.HDF5Methods       as hm 
import PyCSPadImage.GlobalGraphics    as gg
#------------------------------

def test_CSPAD2x2PixCoords() :
    """Test demonstration of how to work with CSPAD2x2CalibPars and CSPAD2x2PixCoords modules
    """    
    #======= Define input parameters
    Ndet = 5
    run  = 180
    path = '/reg/d/psdm/mec/mec73313/calib/CsPad2x2::CalibV1/MecTargetChamber.0:Cspad2x2.%1d/' % Ndet
    #path = '/reg/neh/home1/dubrovin/LCLS/CSPad2x2Alignment/calib-cspad2x2-0%1d-2013-02-13/' % Ndet
    #fname  = '/reg/d/psdm/mec/mec73313/hdf5/mec73313-r%04d.h5' % run
    fname  = '/reg/neh/home1/dubrovin/LCLS/HDF5Analysis-v01/PyCSPadImage/src/mec73313-r%04d.h5' % run
    dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad2x2::ElementV1/MecTargetChamber.0:Cspad2x2.%1d/data' % Ndet
    list_of_clib_types = ['center', 'tilt', 'pedestals']

    #======= Get calibration object
    calib = calpars.CSPAD2x2CalibPars(path, run, list_of_clib_types)

    #======= Get CSPAD2x2 pixel coordinate arrays, shaped as (2, 185, 388)
    coord = pixcoor.CSPAD2x2PixCoords(calib)
    X,Y = coord.get_cspad2x2_pix_coordinate_arrays_pix()

    #======= Get CSPAD2x2 pedestals array, shaped as (185, 388, 2)
    peds_arr = calib.getCalibPars('pedestals', run)

    #======= Get data array from hdf5 dataset, shaped as (185, 388, 2)
    data_arr = hm.getDataSetForOneEvent(fname, dsname, event=0) - peds_arr
    
    #======= Convert shape from (185, 388, 2) to (2, 185, 388)    
    ord_arr  = calpars.data2x2ToTwo2x1(data_arr)

    #======= Compose and plot CSPAD2x2 image from coordinate and intensity arrays
    img2d = gg.getImageFromIndexArrays(X,Y,ord_arr)

    #======= Print for test purpose 
    calib.printCalibParsStatus()
    #print 'pedestals:\n', calib.getCalibPars('pedestals')
    print 'center:\n',    calib.getCalibPars('center')
    print 'tilt:\n',      calib.getCalibPars('tilt')
    print 'peds_arr.shape:', peds_arr.shape  # = (185, 388, 2)  
    print 'Get data array from file: ' + fname
    print 'data_arr.shape:', data_arr.shape
    print 'ord_arr.shape:', ord_arr.shape
    print 'img2d.shape:', img2d.shape

    #======= Plot image and spectrum
    my_range = (-10,40) # None
    gg.plotImageLarge(img2d, amp_range=my_range)        
    gg.plotSpectrum(img2d, amp_range=my_range)
    gg.show()

#------------------------------

def test_cspad2x2_image_with_data() :
    """    
    Example showing how to get calibration, coordinate objects, get pedestals and data arrays and plot image.
    1. get calibration and pixel coordinate objects
    2. get pedestals from calibration
    3. get data for one event from hdf5 file and subtract pedestal array
    4. convert data to image using coord.get_cspad2x2_image(data)
    5. plot image
    """    
    path  = '/reg/d/psdm/mec/mec73313/calib/CsPad2x2::CalibV1/MecTargetChamber.0:Cspad2x2.5/'
    run   = 180
    calib = calpars.CSPAD2x2CalibPars(path, run)
    coord = pixcoor.CSPAD2x2PixCoords(calib)
    peds  = calib.getCalibPars('pedestals', run)

    fname = '/reg/neh/home1/dubrovin/LCLS/HDF5Analysis-v01/PyCSPadImage/src/mec73313-r%04d.h5' % run
    dsname= '/Configure:0000/Run:0000/CalibCycle:0000/CsPad2x2::ElementV1/MecTargetChamber.0:Cspad2x2.5/data'
    data  = hm.getDataSetForOneEvent(fname, dsname, event=0) - peds
    print 'data.shape =', data.shape

    img2d = coord.get_cspad2x2_image(data)
    print 'img2d.shape =', img2d.shape
 
    gg.plotImageLarge(img2d, amp_range=(-10, 40), figsize=(12,11))
    gg.show()

#------------------------------

def test_cspad2x2_image() :
    coord = pixcoor.test_instantiation_2()

#------------------------------
 
if __name__ == "__main__" :
    if len(sys.argv)==1   : print 'Use command: python', sys.argv[0], '<test-number=0-5, 11-16>'
    elif sys.argv[1]=='0' : pixcoor.test_0()
    elif sys.argv[1]=='1' : pixcoor.test_1()
    elif sys.argv[1]=='2' : pixcoor.test_2()
    elif sys.argv[1]=='3' : pixcoor.test_3()
    elif sys.argv[1]=='4' : test_cspad2x2_image_with_data()
    elif sys.argv[1]=='5' : test_CSPAD2x2PixCoords()

    elif sys.argv[1]=='11': calpars.main_test()
    elif sys.argv[1]=='12': calpars.test_reshaping_arrs_for_cspad2x2()
    elif sys.argv[1]=='13': cpd.main_test()
    elif sys.argv[1]=='14': pixcoor2x1.test_2x1_xy_maps()
    elif sys.argv[1]=='15': pixcoor2x1.test_2x1_img()
    elif sys.argv[1]=='16': pixcoor2x1.test_2x1_img_easy()()
    
    else : print 'Non-expected arguments: sys.argv=', sys.argv

    sys.exit ( 'End of test.' )

#------------------------------

