#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module HDF5Methods...
#
#------------------------------------------------------------------------

"""This module provides access to the calibration parameters

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
import sys
import os
import h5py
import numpy as np

#import ConfigParameters as cp
import CSPadConfigPars  as ccp
import GlobalMethods    as gm

#---------------------
#  Class definition --
#---------------------

class HDF5File(object) :
    """This class contains a few methods to manipulate with hdf5 files"""

#---------------------

    def __init__ (self, fname=None) :
        #print """HDF5File: Initialization"""
        self.dset  = None   
        self.fname = fname   
        if fname != None : h5file = self.open_hdf5_file(fname)
        else             : self.h5file = None

#---------------------

    def open_hdf5_file(self, fname) :

        #print '=== Open HDF5 file: ' + fname
        if not os.path.exists(fname) :
            print  'ERROR: THE SPECIFIED FILE "', fname, '" DOES NOT EXIST.'
            return None
            #sys.exit ( 'Exit on ERROR' )

        try :
            self.h5file = h5py.File(fname,'r') # open read-only
            self.fname = fname
            print  'Open file', fname
            return self.h5file

        except IOError:
            print 'IOError: CAN NOT OPEN FILE:', fname
            return None
            #sys.exit ( 'Exit on ERROR' )

#---------------------

    def close_hdf5_file(self) :
        self.h5file.close()
        #print '=== Close HDF5 file ==='

#---------------------

    def get_dataset_from_hdf5_file(self,dsname) :
        #print 'From hdf5 file get dataset :', dsname
        try :
            self.dset = self.h5file[str(dsname)]
            return self.dset
        except KeyError:
            #print 'ERROR: DATASET %s \nDOES NOT EXIST IN HDF5 file %s' % (dsname, self.fname)
            return None
            #sys.exit ( 'Exit on ERROR' )

#---------------------

    def get_cspad_config_dsname( self, data_dsname ) :
        """Find the CSPAD configuration dataset name in hdf5 file."""

        grpname = '/Configure:0000'
        pattern = 'CsPad::ConfigV'
        suffix  = gm.get_item_second_to_last_name(data_dsname) + '/config' 

        #print 'get_cspad_config_dsname(): loop over group content:'
        grp = self.get_dataset_from_hdf5_file(grpname)
        for key,val in dict(grp).iteritems() :
            #print '  ', key, val
            if key.find(pattern)==0 :
                #print '    ', val.name
                dsname = val.name  + '/' + suffix
                #print 'get_cspad_config_dsname(): found configuration dsname in hdf5 file:', dsname
                return dsname

        return None

#---------------------


#---------------------
#---------------------
# Out of class methods
#---------------------
#---------------------

def print_hdf5_item_structure(g, offset='    ') :
    """Prints the input file/group/dataset (g) name and begin iterations on its content"""
    if   isinstance(g,h5py.File) :
        print g.file, '(File)', g.name

    elif isinstance(g,h5py.Dataset) :
        print '(Dataset)', g.name, '    len =', g.shape #, g.dtype

    elif isinstance(g,h5py.Group) :
        print '(Group)', g.name

    else :
        print 'WORNING: UNKNOWN ITEM IN HDF5 FILE', g.name
        sys.exit ( "EXECUTION IS TERMINATED" )

    if isinstance(g, h5py.File) or isinstance(g, h5py.Group) :
        for key,val in dict(g).iteritems() :
            subg = val
            print offset, key, #,"   ", subg.name #, val, subg.len(), type(subg),
            print_hdf5_item_structure(subg, offset + '    ')

#---------------------

def get_cspad_name_and_data_type_from_dataset_name( dsname ) :
    """Returns the detector name and data type as a parts of the dataset name
    """
    cspad_name      = ''
    cspad_data_type = ''
    if gm.CSpadIsInTheName(dsname) :
        name1, name2, name3 = gm.get_item_last_three_names(dsname)
        cspad_name      = name2 # something like CxiDs1.0:Cspad.0
        cspad_data_type = name3 # something like CsPad::ElementV2
        #print 'CSPad is found in dataset', dsname
        #print 'CSPad name and data type =', cs.confcspad.cspad_name, cs.confcspad.cspad_data_type
        #return name2
    else :
        print 'ERROR in get_cspad_name_and_data_type_from_dataset_name( dsname ) for\ndsname =', dsname,\
              '\nReturn empty the CSPad name and data_type !!!'
    return cspad_name, cspad_data_type

#---------------------

def get_root_calib_dir_from_hdf5_path( path ) :
    """Returns the root calib directory, i.e. /reg/d/psdm/CXI/cxi80410/calib/.
       Search for the "/hdf5" in the path or full file name, split it and add the "/calib/".
    """
    ind1 = path.find('/hdf5')
    if ind1==-1 :
        print 'ERROR in get_root_calib_dir_from_hdf5_path(): The input path:', path, \
        'does not contain the expected "/hdf" subdir !!! Return empty root path.'
        return ''
        
    return path[0:ind1] + '/calib/'

#---------------------

def get_cspad_calib_dir( path, cspad_name='CxiDs1.0:Cspad.0' ) :
    """Constructs the calibration directory name from the path to *.h5 file and the detector name.
       The lates version of calibration is picked up automatically.
    """
    root_calib_dir = get_root_calib_dir_from_hdf5_path( path )

    dirList=gm.getListOfFilesInDir(root_calib_dir)
    calib_version = '' 
    for name in dirList: # loop
        if name[0:5] == 'CsPad' :
            calib_version = name # something like CsPad::CalibV1

    return root_calib_dir + calib_version + '/' + cspad_name     

#---------------------

def get_run_number_from_hdf5_file_name( fname ) :
    """Returns the run number from the *-r<run_number>.h5 file name
    """
    ind1 = fname.find('-r')
    ind2 = fname.find('.h5')
    if ind1==-1 or ind2==-1 or ind2<=ind1 :
        print 'ERROR in get_run_number_from_hdf5_file_name(): The input file name:', fname, \
        'is not an expected *-r<run_number>.h5 file name (i.e. cxi80410-r0742.h5) !!! Return the run_num=0'
        return 0
    #print 'File name:', fname, 'Run number:', cs.confcspad.run_num
    return int(fname[ind1+2:ind2])

#---------------------

def getDataSetForOneEvent( fname  = '/reg/d/psdm/CXI/cxi35711/hdf5/cxi35711-r0009.h5',
                           dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/CxiDs1.0:Cspad.0/data',
                           event  = 0 ) :

    #print 'fname:', fname
    #print 'dsname:', dsname

    hdf5file = hdf5mets.open_hdf5_file(fname)
    dataset  = hdf5mets.get_dataset_from_hdf5_file(dsname)
    evdata   = dataset[event]
    hdf5mets.close_hdf5_file()
    return evdata

#---------------------

def getOneCSPadEventForTest( fname  = '/reg/d/psdm/CXI/cxi35711/hdf5/cxi35711-r0009.h5',
                             dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/CxiDs1.0:Cspad.0/data',
                             event  = 0 ) :

    #ccp.cspadconfig.setCSPadConfiguration( fname, dsname, event )
    #ccp.cspadconfig.printCSPadConfigPars()

    #file    = h5py.File(fname, 'r')
    #dataset = file[dsname]
    #evdata  = dataset[event]
    #file.close()
    #print 'fname:', fname
    #print 'dsname:', dsname
    
    h5file = hdf5mets.open_hdf5_file(fname)
    ccp.cspadconfig.setCSPadConfigurationFromOpenFile( hdf5mets, dsname, event )
    dataset  = hdf5mets.get_dataset_from_hdf5_file(dsname)
    evdata   = dataset[event]
    hdf5mets.close_hdf5_file()
    return evdata

#---------------------

def getAverageCSPadEvent( fname   = '/reg/d/psdm/CXI/cxi35711/hdf5/cxi35711-r0009.h5',
                          dsname  = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/CxiDs1.0:Cspad.0/data',
                          event1  = 1,
                          nevents = 10) :

    print 'Average over', nevents, 'events, starting from', event1

    hdf5file = hdf5mets.open_hdf5_file(fname)
    ccp.cspadconfig.setCSPadConfigurationFromOpenFile( hdf5mets, dsname, event1 )
    dataset  = hdf5mets.get_dataset_from_hdf5_file(dsname)

    evdata = np.zeros(dataset[event1].shape, dtype=np.float32)
    for evt in range(event1, event1+nevents) :
        evdata += dataset[evt]
        if evt%10 == 0 :
            print 'add event', evt
            #print evdata[1,:]

    evdata /= nevents

    hdf5mets.close_hdf5_file()
    return evdata

#---------------------

hdf5mets = HDF5File()

#----------------------------------------------
#----------------------------------------------
#----------------------------------------------
#----------------------------------------------

def main_test1() :

    event   = 1
    fname   = '/reg/d/psdm/CXI/cxi35711/hdf5/cxi35711-r0009.h5'
    dsname  = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/CxiDs1.0:Cspad.0/data'

    h5file = hdf5mets.HDF5File(fname)

    grp = h5file.get_dataset_from_hdf5_file('/')    
    print_hdf5_item_structure(grp)

    arr = h5file.get_dataset_from_hdf5_file(dsname)
    print 'arr[event]=\n', arr[event]

    h5file.close_hdf5_file()

#---------------------

def main_test2() :
    
    fname = 'cxi80410-r019742.h5'
    print '\nTest get_run_number_from_hdf5_file_name(', fname, '), \nrun=',\
          get_run_number_from_hdf5_file_name(fname)

    dsname = '/Configure:0000/Run:0000/CalibCycle:0000/CsPad::ElementV2/CxiDs1.0:Cspad.0/data'
    print '\nTest get_cspad_name_and_data_type_from_dataset_name(', dsname,') \nCSPad (name, data_type) =',\
          get_cspad_name_and_data_type_from_dataset_name(dsname)

    path_and_fname = '/reg/d/psdm/CXI/cxi35711/hdf5/cxi35711-r0009.h5'
    print '\nTest get_root_calib_dir_from_hdf5_path(', path_and_fname,') \nroot_calib_dir =',\
          get_root_calib_dir_from_hdf5_path( path_and_fname )

    print '\nTest get_cspad_calib_dir(', path_and_fname,') \ncspad_calib_dir =',\
          get_cspad_calib_dir( path_and_fname, cspad_name='CxiDs1.0:Cspad.0' )

#---------------------

if __name__ == "__main__" :

    #main_test1()
    main_test2()
    sys.exit ( 'End of test.' )

#----------------------------------------------
