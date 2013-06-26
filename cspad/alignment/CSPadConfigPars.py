#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module CSPadConfigPars...
#
#------------------------------------------------------------------------

"""This module provides access to the CSPad configuration parameters

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
    """This class provides access to the CSPad configuration parameters:
    1. Sets default configuration parameters;
    2. Sets actual coniguration parameters using the hdf5 file, so it needs in
         * hdf5 file name:       fname or cp.confpars....
         * dsname for CSPad data
         * current event number: event or cp.confpars.eventCurrent
    Most important parameters, which may have changed from run-to-run are:
    indPairsInQuads,
    quadNumsInEvent
    """
    nquads      =   4    # Total number of quads in cspad
    nsects      =   8    # Total number of sections in quad
    wid2x1      = 185    # Number of rows in 2x1 at rotation 0
    len2x1      = 388    # Number of cols in 2x1 at rotation 0
    gapIn2x1    =   3    # Gap in 2x1 between two halfs

    pixSize     = 109.92 # Pixel size in um (micrometer)
    pixSizeWide = 274.80 # Wide pixel size in um (micrometer)


#---------------------

    def __init__ (self) :
        """Define default parameters"""

        #self.indPairsInQuads = np.arange(32)
        #self.indPairsInQuads.shape = (4,8)
        self.indPairsInQuads = np.array( 
            [ [ 0,   1,   2,   3,   4,   5,   6,   7],
              [ 8,   9,  10,  11,  12,  13,  14,  15],
              [16,  17,  18,  19,  20,  21,  22,  23],
              [24,  25,  26,  27,  28,  29,  30,  31] ] )

        self.quadNumsInEvent = [0, 1, 2, 3] # Numeration of quads in the event record

        self.pairInQaudOriInd = np.array(
            [ [   3,   3,   2,   2,   1,   1,   2,   2],
              [   3,   3,   2,   2,   1,   1,   2,   2],
              [   3,   3,   2,   2,   1,   1,   2,   2],
              [   3,   3,   2,   2,   1,   1,   2,   2] ])

        self.quadInDetOriInd = [2, 1, 0, 3]

        self.margin    = 4    # For tilt spare space near frame

        self.cspadQuad = 0    # Defauld quad number
        self.cspadPair = 0    # Defauld pair/section number

        self.quadDimX = 850   # Quad image X dimension 
        self.quadDimY = 850   # Quad image Y dimension 

        self.detDimX  = 1765  # Quad image X dimension 
        self.detDimY  = 1765  # Quad image Y dimension 

        self.isCSPad2x2 = False # distinguish the CSPad2x2 (CSPad2x2Element) from regular CSPad (CSPadElement)

       #self.dsnameCSpadV2    = "/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/XppGon.0:Cspad.0/data"
       #self.dsnameCSpadV3Conf= "/Configure:0000/CsPad::ConfigV3/"                        #CxiDs1.0:Cspad.0/config - is added auto
        self.dsnameCSpadVXConf= "/Configure:0000/CsPad::ConfigV"

       #self.fileNameWithAlreadySetCSPadConfiguration = 'File name for CSPad configuration is not set yet...'

#---------------------

    def getQuadNumsInEvent( self, dsname, event=0 ):
        """For each event"""
        el_dsname  = gm.get_item_path_to_last_name(dsname) + '/element'
        try:
            ds_element = self.h5file[el_dsname]
        except KeyError: 
            print 80*'!'
            print 'WARNING: The CSPad configuration for "element" dataset is not found. Default will be used'
            print 80*'!'
            return self.quadNumsInEvent

        try:
            return ds_element[event]['quad']
        except KeyError: 
            print 80*'!'
            print 'WARNING: The CSPad configuration for ds_element[event][quad] is not found. Default will be used'
            print 80*'!'
            return self.quadNumsInEvent

#---------------------

    def getIndPairsInQuads( self, dsname ):
        """For each run"""
        item_second_to_last_name = gm.get_item_second_to_last_name(dsname)
        #Find the configuration dataset name for CSPad
        self.dsConf = None
        for vers in range(100) :
            self.cspad_config_ds_name = self.dsnameCSpadVXConf + str(vers) + '/' + item_second_to_last_name + '/config' 

            #print 'Try to get CSPad configuration from the dataset:\n',self.cspad_config_ds_name 

            try:
                self.dsConf = self.h5file[self.cspad_config_ds_name]      # t=0.01us
                break
            except KeyError:
                #print 'Try another configuration dataset name'
                continue

        if self.dsConf == None :
            print 80*'!'
            print 'WARNING: The CSPad configuration for "config" dataset is not found. Default will be used'
            print 80*'!'
            return self.indPairsInQuads
        else :
           #print 'Configuration dataset name:', self.cspad_config_ds_name
            return self.dsConf['sections'] # For V2 it is dsConf.value[13], for V3 it is 15th  ...

#---------------------

    def setCSPadConfiguration( self, fname, dsname, event=0 ):
        """Takes the CSPad configuration parameters from hdf5 file."""
        if gm.CSpad2x2ElementIsInTheName(dsname) :
            print 'getCSpadConfiguration(...): This is a CSpad2x2Element. Special configuration is not required'
            self.isCSPad2x2 = True
            return

        self.h5file = hm.hdf5mets.open_hdf5_file(fname) 

        self.quadNumsInEvent = self.getQuadNumsInEvent( dsname, event )
        self.indPairsInQuads = self.getIndPairsInQuads( dsname )

        #if fname != self.fileNameWithAlreadySetCSPadConfiguration :
        # Once per file:
            #self.fileNameWithAlreadySetCSPadConfiguration = fname            
            #print "Indexes of pairs in quads =\n", self.indPairsInQuads 

        hm.hdf5mets.close_hdf5_file()
        #self.printCSPadConfigPars()

#---------------------

    def setCSPadConfigurationFromOpenFile( self, h5file, dsname, event=0 ):
        """Takes the CSPad configuration parameters from open hdf5 file."""
        if gm.CSpad2x2ElementIsInTheName(dsname) :
            print 'getCSpadConfiguration(...): This is a CSpad2x2Element. Special configuration is not required'
            self.isCSPad2x2 = True
            return

        self.h5file = h5file
        self.quadNumsInEvent = self.getQuadNumsInEvent( dsname, event )
        self.indPairsInQuads = self.getIndPairsInQuads( dsname )
        #self.printCSPadConfigPars()

#---------------------

    def printCSPadConfigPars(self) :
        print 50*'-'
        print 'printCSPadConfigPars():'
        print '\nindPairsInQuads  =\n', self.indPairsInQuads
        print '\nquadNumsInEvent  =',   self.quadNumsInEvent
        print 50*'-'

#---------------------

cspadconfig = CSPadConfigPars() 

#----------------------------------------------

def main() :
    """This is an example of how to use this class"""
    #fname  = '/reg/d/psdm/CXI/cxi35711/hdf5/cxi35711-r0009.h5'
    #dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/CxiDs1.0:Cspad.0/data'
    #event  = 1

    fname  = '/reg/d/psdm/CXI/cxi37411/hdf5/cxi37411-r0039.h5'
    dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/CxiDsd.0:Cspad.0/data'
    event  = 1

    print 'Default CSPad configuration pars:'
    cspadconfig.printCSPadConfigPars()

    print '\nCSPad configuration pars: for fname, dsname, event =\n', fname, '\n', dsname, '\n', event
    cspadconfig.setCSPadConfiguration( fname, dsname, event ) # This will set current CSPad configuration
    cspadconfig.printCSPadConfigPars()


if __name__ == "__main__" :
    main()
    sys.exit ( 'End of test.' )

#----------------------------------------------
