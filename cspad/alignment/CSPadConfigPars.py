#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module CSPadConfigPars...
#
#------------------------------------------------------------------------

"""This module provides access to the CSPAD configuration parameters

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
#import os
import sys

import numpy            as np
import GlobalMethods    as gm
import HDF5Methods      as hm
#import ConfigParameters as cp

#---------------------
#  Class definition --
#---------------------

class CSPadConfigPars(object) :
    """This class provides access to the CSPAD configuration parameters.

    1. Sets a bunch of default configuration parameters,
       loads current coniguration parameters from hdf5 file or from external parameters.
    2. Provides access to current coniguration parameters
    3. Contains conversion methods for arrays between raw data and entire cspad.
    
    Interface
    =========
    1.  Instatiation
    1.1 Default constructor sets default parameters for indPairsInQuads & quadNumsInEvent:

        config = CSPadConfigPars()

    1.2 Initialization of configuration parameters using hdf5 file, for example:
        fname  = '/reg/d/psdm/xpp/xpp66213/hdf5/xpp66213-r0150.h5'
        dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/XppGon.0:Cspad.0/data'

        config.setCSPadConfiguration( fname, dsname, event=0 ):

    1.3 Initialization of configuration parameters using external arrays, for example:        
        indPairs = np.arange(32)
        indPairs.shape = (4,8)
        quadNums = np.arange(4)

        config.setCSPadConfigArrays( indPairsInQuads=indPairs, quadNumsInEvent=quadNums )

    2.  Access methods:

    2.1 Access to indPairsInQuads and quadNumsInEvent:
        quadNums = config.getQuadNumsInEvent()
        indPairs = config.getIndPairsInQuads()
        config.printCSPadConfigPars()

    2.2 Access to static class parameters:    
        import CSPadConfigPars as ccp
        my_wid2x1 = ccp.CSPadConfigPars().wid2x1
        etc...

    3.  Conversions between entire (4,8,185,388) and shaped as data (N<32,185,388) cspad pixel array shapes:

    3.1 Conversion of the entire cspad pixel array arr_entire_cspad with shape (4,8,185,388)
        in to the arr_raw_data, shaped as data (N<32,185,388):
        arr_raw_data = config.getCSPadPixArrayShapedAsData(arr_entire_cspad)

    3.2 Conversion of the cspad pixel array arr_raw_data shaped as data (N<32,185,388)
        in to the entire cspad pixel array arr_entire_cspad with shape (4,8,185,388):
        arr_entire_cspad = getCSPadPixArrayFromArrayShapedAsData(arr_raw_data)

    4.  Tests
        To test CSPadConfigPars from release directory use command:
        python PyCSPadImage/src/CSPadConfigPars.py <test-number>
        where <test-number> stands for  0, 1, 2, or 3
    """

    nquads      =   4    # Total number of quads in cspad
    nsects      =   8    # Total number of sections in quad
    wid2x1      = 185    # Number of rows in 2x1 at rotation 0
    len2x1      = 388    # Number of cols in 2x1 at rotation 0
    gapIn2x1    =   3    # Gap in 2x1 between two halfs

    pixSize     = 109.92 # Pixel size in um (micrometer)
    pixSizeWide = 274.80 # Wide pixel size in um (micrometer)

    pairInQaudOriInd = np.array(
            [ [   3,   3,   2,   2,   1,   1,   2,   2],
              [   3,   3,   2,   2,   1,   1,   2,   2],
              [   3,   3,   2,   2,   1,   1,   2,   2],
              [   3,   3,   2,   2,   1,   1,   2,   2] ])

    quadInDetOriInd = [2, 1, 0, 3]

    margin    = 4    # For tilt spare space near frame

    cspadQuad = 0    # Defauld quad number
    cspadPair = 0    # Defauld pair/section number

    quadDimX = 850   # Quad image X dimension 
    quadDimY = 850   # Quad image Y dimension 

    detDimX  = 1765  # Quad image X dimension 
    detDimY  = 1765  # Quad image Y dimension 

    isCSPad2x2 = False # distinguish the CSPad2x2 (CSPad2x2Element) from regular CSPad (CSPadElement)

    #data_dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/XppGon.0:Cspad.0/data'
    #conf_dsname = '/Configure:0000/CsPad::ConfigV4                          /XppGon.0:Cspad.0/config'

    #indPairsInQuads_def = np.arange(32)
    #indPairsInQuads_def.shape = (4,8)
    indPairsInQuads_def = np.array( 
            [ [ 0,   1,   2,   3,   4,   5,   6,   7],
              [ 8,   9,  10,  11,  12,  13,  14,  15],
              [16,  17,  18,  19,  20,  21,  22,  23],
              [24,  25,  26,  27,  28,  29,  30,  31] ] )

    quadNumsInEvent_def = [0, 1, 2, 3] # Numeration of quads in the event record

#---------------------

    def __init__ (self) :
        """Sets default parameters"""

        self.setCSPadConfigArrays() # set default configuration

#---------------------

    def setCSPadConfigArrays( self, indPairsInQuads=None, quadNumsInEvent=None ):
        """Sets the CSPAD configuration from input parameters."""

        self.ds_element = None

        if indPairsInQuads != None : self.indPairsInQuads = indPairsInQuads
        else                       : self.indPairsInQuads = self.indPairsInQuads_def
        
        if indPairsInQuads != None : self.quadNumsInEvent = quadNumsInEvent
        else                       : self.quadNumsInEvent = self.quadNumsInEvent_def

#---------------------

    def _setQuadNumsInEvent( self, h5file, dsname, event=0 ):
        """Sets the self.quadNumsInEvent from h5file using dsname and optional event number."""

        el_dsname  = gm.get_item_path_to_last_name(dsname) + '/element'
        #self.ds_element = h5file.get_dataset_from_hdf5_file(el_dsname)[:]['quad']
        self.ds_element = h5file.get_dataset_from_hdf5_file(el_dsname)
        #print 'self.ds_element = ', self.ds_element
        self.quadNumsInEvent = self.getQuadNumsInEvent( event )

#---------------------

    def getQuadNumsInEvent( self, event=0 ):
        """Gets the uadNumsInEvent from self.ds_element dataset (if available) for given event number."""

        if event==0 or self.ds_element == None : return self.quadNumsInEvent

        try:
            #return self.ds_element[event]
            return self.ds_element[event]['quad']
        except KeyError: 
            print 80*'!'
            print 'WARNING: The CSPad configuration for ds_element[event][quad] is not found. Default will be used'
            print 80*'!'
            return self.quadNumsInEvent_def

#---------------------

    def getIndPairsInQuads( self ):
        """Returns current np.array for indPairsInQuads."""

        return self.indPairsInQuads
        
#---------------------

    def _setIndPairsInQuads( self, h5file, dsname ):
        """Sets the self.indPairsInQuads from h5file using dsname."""

        config_dsname = h5file.get_cspad_config_dsname(dsname)        
        config_ds     = h5file.get_dataset_from_hdf5_file(config_dsname)

        try:
            self.indPairsInQuads = config_ds['sections'] # For V2 it is config_ds.value[13], for V3 it is 15th  ...
        except : # NoneType KeyError: 
            print 80*'!'
            print 'WARNING: The CSPAD configuration for "config" dataset is not found. Default will be used'
            print 80*'!'
            self.indPairsInQuads = self.indPairsInQuads_def

#---------------------

    def _isCSPad2x2Dataset( self, dsname ):
        """Check if the dsname is for CsPad2x2."""
        if gm.CSpad2x2ElementIsInTheName(dsname) :
            print 'getCSpadConfiguration(...): This is a CSpad2x2Element. Special configuration is not required'
            self.isCSPad2x2 = True
        else:
            self.isCSPad2x2 = False
        return self.isCSPad2x2

#---------------------

    def setCSPadConfigurationFromOpenFile( self, h5file, dsname, event=0 ):
        """Sets the CSPad configuration parameters from open hdf5 file."""

        if self._isCSPad2x2Dataset(dsname) : return

        self._setQuadNumsInEvent( h5file, dsname, event )
        self._setIndPairsInQuads( h5file, dsname )

#---------------------

    def setCSPadConfiguration( self, fname, dsname, event=0 ):
        """Sets the CSPAD configuration parameters from hdf5 file."""

        h5file = hm.HDF5File(fname) 
        self.setCSPadConfigurationFromOpenFile( h5file, dsname, event )
        h5file.close_hdf5_file()

#---------------------

    def getCSPadPixArrayShapedAsData(self, arr_in):
        """Converts the entire cspad pixel array arr_in with shape (4,8,185,388)
           in to the arr_out, shaped as data (nsects_in_data,185,388).
        """

        if arr_in.shape != (4,8,185,388) :
            try :
                arr_in.shape = (4,8,185,388)
            except :
                print 'ERROR in getCSPadPixArrayShapedAsData(): Input array shape=', arr_out.shape,\
                      'can not be reshaped to (4,8,185,388)'
                return arr_in

        nsects_in_data = max(self.indPairsInQuads.flatten()) + 1
        arr_out = np.zeros((nsects_in_data,185,388), dtype=arr_in.dtype) # dtype=np.float32)
          
        for iq in range(len(self.quadNumsInEvent)) :
            quad = int(self.quadNumsInEvent[iq]) # uint8 -> int
            for sect in range(8): # loop over ind = 0,1,2,...,7
                ind_sect_in_data = self.indPairsInQuads[quad][sect]
                if ind_sect_in_data == -1 : continue
                arr_out[ind_sect_in_data][:] = arr_in[quad][sect][:]

        print 'getCSPadPixArrayShapedAsData(): Created arr_out.shape=', arr_out.shape
        return np.array(arr_out) # creates a copy of the pixel array
 
#---------------------

    def getCSPadPixArrayFromArrayShapedAsData(self, arr_in):
        """Converts the cspad pixel array shaped as data (nsects_in_data,185,388).
           in to the entire cspad pixel array arr_out with shape (4,8,185,388)
        """

        nsects_in_data = max(self.indPairsInQuads.flatten()) + 1
        arr_out = np.zeros((4,8,185,388), dtype=arr_in.dtype) # dtype=np.float32)
          
        for iq in range(len(self.quadNumsInEvent)) :
            quad = int(self.quadNumsInEvent[iq]) # uint8 -> int
            for sect in range(8): # loop over ind = 0,1,2,...,7
                ind_sect_in_data = self.indPairsInQuads[quad][sect]
                if ind_sect_in_data == -1 : continue

                arr_out[quad][sect][:] = arr_in[ind_sect_in_data][:]

        print 'getCSPadPixArrayFromArrayShapedAsData(): Created arr_out.shape=', arr_out.shape
        return np.array(arr_out) # creates a copy of the pixel array
 
#---------------------

    def printCSPadConfigPars(self) :
        """Prints the CSPad current configuration parameters."""
        print 50*'-'
        print 'printCSPadConfigPars():'
        print '\nquadNumsInEvent()  =',   self.getQuadNumsInEvent()
        print '\nindPairsInQuads()  =\n', self.getIndPairsInQuads()
        print 50*'-'

#---------------------

cspadconfig = CSPadConfigPars() 

#----------------------------------------------

def get_test_fname_dsname_event() :
    """This is an example of how to use this class"""
    #fname  = '/reg/d/psdm/CXI/cxi35711/hdf5/cxi35711-r0009.h5'
    #fname  = '/reg/d/psdm/CXI/cxi37411/hdf5/cxi37411-r0039.h5'
    fname  = '/reg/d/psdm/xpp/xpp66213/hdf5/xpp66213-r0150.h5'

    #dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/CxiDsd.0:Cspad.0/data'
    #dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/CxiDsd.0:Cspad.0/data'
    dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/XppGon.0:Cspad.0/data'

    event  = 1

    return fname, dsname, event

#----------------------------------------------

from time import time # for test purpose only

def test_object(config) :
    """Test the config object using overloading of parameters from hdf5 file."""

    fname, dsname, event  = get_test_fname_dsname_event()

    print 'Default CSPad configuration pars:'
    config.printCSPadConfigPars()

    print '\nCSPad configuration pars: for fname, dsname, event =\n', fname, '\n', dsname, '\n', event
    t0_sec = time()
    config.setCSPadConfiguration( fname, dsname, event ) # This will set current CSPad configuration
    print 'Consumed time to get CSPAD configuration parameters from hdf5 file (sec) =', time()-t0_sec
    config.printCSPadConfigPars()

#----------------------------------------------

def test_object_with_external_pars() :
    """Test of CSPadConfigPars object with initialization using external parameters"""

    config = CSPadConfigPars() # instatiate object

    indPairs = np.array( 
        [ [ 0,   1,   2,   3,   4,   5,   6,   7],
          [ 8,   9,  10,  11,  12,  13,  14,  15],
          [16,  17,  18,  19,  20,  21,  22,  23],
          [24,  -1,  25,  26,  27,  -1,  29,  30] ] )
    quadNums = [2, 3, 0, 1]
    config.setCSPadConfigArrays( indPairsInQuads=indPairs, quadNumsInEvent=quadNums )
    config.printCSPadConfigPars()

    print '\nTest of conversion methods'
    print '\n1. Conversion from entire cspad array to shaped as raw data'

    arr_entire_cspad = np.arange(4*8*185*388)
    arr_entire_cspad.shape = (4,8,185,388)
    t1_sec = time()
    arr_raw_data = config.getCSPadPixArrayShapedAsData(arr_entire_cspad)
    print 'getCSPadPixArrayShapedAsData transformation time (sec) =', time()-t1_sec

    #print 'arr_raw_data:\n', arr_raw_data
    print 'arr_raw_data.shape:', arr_raw_data.shape
    
    print '\n2. Conversion from raw data to entire cspad array'
    t2_sec = time()
    arr_entire_cspad_2 = config.getCSPadPixArrayFromArrayShapedAsData(arr_raw_data)
    print 'getCSPadPixArrayFromArrayShapedAsData transformation time (sec) =', time()-t2_sec

    #print 'arr_entire_cspad_2:\n',     arr_entire_cspad_2
    print 'arr_entire_cspad_2.shape:', arr_entire_cspad_2.shape

    print '\n3. Check if arrays are equals after two transformations'
    indPairsInQuads = config.getIndPairsInQuads()
    print 'indPairsInQuads:\n', indPairsInQuads
    arraysAreEqual = True
    for q in range(4) :
        for s in range(8) :
            if indPairsInQuads[q,s] != -1 :
                if not np.equal(arr_entire_cspad[q,s,:],arr_entire_cspad_2[q,s,:]).all() : arraysAreEqual = False

    if arraysAreEqual : print 'Arrays are equal in expected parts after two transformations'
    else              : print 'Arrays are NOT equal after two transformations'

#----------------------------------------------

if __name__ == "__main__" :

    print 65*'='

    if len(sys.argv)==1 :
        print 'Use command: python', sys.argv[0], '<test-number=0-3>'

    elif sys.argv[1]=='0' :
        print 'Test singleton of CSPadConfigPars'
        config = cspadconfig # This is a singleton name
        test_object(config)

    elif sys.argv[1]=='1' :
        print 'Test default CSPadConfigPars object'
        config = CSPadConfigPars() # instatiate object
        test_object(config)

    elif sys.argv[1]=='2' :
        print 'Test of access to static CSPadConfigPars parameters'
        print 'my_wid2x1:', CSPadConfigPars().wid2x1
        print 'my_len2x1:', CSPadConfigPars().len2x1

    elif sys.argv[1]=='3' :
        print 'Test of CSPadConfigPars object with initialization using external parameters'
        test_object_with_external_pars()

    else :
        print 'Non-expected arguments: sys.argv=', sys.argv

    sys.exit ( 'End of test.' )

#----------------------------------------------
