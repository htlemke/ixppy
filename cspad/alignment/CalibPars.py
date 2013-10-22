#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module CalibPars...
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
#from PyQt4 import QtGui, QtCore
import numpy as np
import CalibParsDefault   as cpd
import CSPADCalibParsEvaluated as cpe

#---------------------
#  Class definition --
#---------------------

class CalibPars (object) :
    """This class provides access to the CSPAD calibration parameters.

       Interface
       =========
       Regular instantiation:
           ALL parameters are OPTIONAL NAMED parameters;
           path  = '/reg/d/psdm/xpp/xpptut13/calib/CsPad::CalibV1/XppGon.0:Cspad.0/' 
           run   = 123
           calib = CalibPars(path, list_of_clib_types=['center', 'tilt', 'pedestals'])
           arr_pedestals = calib.getCalibPars('pedestals', run)
 
       Other option for instantiation:
           calib    = CalibPars()
           run      = 123  - is an optional, named
           calibdir = '/reg/d/psdm/CXI/cxi35711/calib'
           group    = 'CsPad::CalibV1'
           source   = 'CxiDs1.0:Cspad.0'
           calib.setCalibPars (run, calibdir, group, source)

       Get array of calibration parameters for specified type and run number:
           type = 'center'
           arr  = calib.getCalibPars (type[,run])
    """

    list_of_clib_types_total = cpd.calibparsdefault.list_of_clib_types
        #[
        # 'center'
        #,'center_corr'  
        #,'marg_gap_shift' 
        #,'offset'
        #,'offset_corr'
        #,'rotation'
        #,'tilt'
        #,'quad_rotation'
        #,'quad_tilt'
        #,'center_global'  
        #,'tilt_global'  
        #,'beam_vector'
        #,'beam_intersect'
        #,'common_mode'
        #,'pedestals'
        #,'filter'
        #,'pixel_status'
        #]

#---------------------

    def __init__ (self, path=None, run=None, list_of_clib_types=None) :
        """Class constructor:
           The path and run need to be defined to instantiate the object and load correct set of parameters.
           If path or run is omitted, default parameters will be used. 
           Run number can be omitted here and passed later in getCalibPars(type, run)
           list_of_clib_types - optional parameter for optimization of time; only types from the list will be loaded.
           If the list_of_clib_types is omitted, all parameters will be loaded and used.
        """
        self.cpars        = {}
        self.cpars_status = {} # 'DEFAULT', 'FROM_FILE', 'EVALUATED'

        self.cpeval = None # Will be defined after input of all parameters

        self.run = run # 10
        self.path_to_calib_types = path 

        self.setListOfCalibTypes(list_of_clib_types)
        
        self.setCalibParsDefault()

        #if path!=None and run!=None :
        #    self.setCalibParsForPath(run, path)

#---------------------

    def setListOfCalibTypes(self, list_of_clib_types) :
        """Defines the internal list of calibration types which will be used.
        """
        if list_of_clib_types == None :
            self.list_of_clib_types = self.list_of_clib_types_total
            return

        self.list_of_clib_types = []
        for type in list_of_clib_types :
            if type in self.list_of_clib_types_total :
                self.list_of_clib_types.append(type)
            else :
                msg = 'WARNING: TYPE ' + type + ' IS UNKNOWN FOR CSPAD2x2' + \
                      '\n KNOWN TYPES:' + str(self.list_of_clib_types_total)
                print msg

#---------------------

    def setCalibParsDefault (self) :
        """Loads default calibration parameters from singleton object.
        """    
        #print 'Set default calibration parameters'
        for type in self.list_of_clib_types :
            self.cpars[type] = self.getCalibParsDefault (type)
            self.cpars_status[type] = 'DEFAULT'

#---------------------

    def getCalibParsDefault (self, type) :
        """Returns the default calibration parameters for specified type.
        """
        return cpd.calibparsdefault.getCalibParsDefault (type)

#---------------------

    def setRun ( self, run=None ) :
        """Resets the run number and loads calibration parameters, if calib files are available.
        """    
        if run!=None  :
            self.run = run
            #print 'Load the calibration parameters for run ', self.run
            #self.setCalibParsForPath()

#---------------------

    def setCalibPars (self,
                      run      = None,   # 1
                      calibdir = None,   # '/reg/d/psdm/CXI/cxi35711/calib',
                      group    = None,   # 'CsPad::CalibV1',
                      source   = None) : # 'CxiDs1.0:Cspad.0') :
        """Set path and run for calibration parameters."""

        #path = os.path.join(calibdir,group,source)
        
        self.setCalibParsForPath (run, calibdir + '/' + group + '/' + source)

#---------------------

    def setCalibParsForPath (self,
                             run  = None,   # 1
                             path = None) : #'/reg/d/psdm/CXI/cxi35711/calib/CsPad::CalibV1/CxiDs1.0:Cspad.0' ) :
        """Set path and run for calibration parameters."""

        self.path_to_calib_types = path
        if run!=None  : self.run = run

        # self.loadAllCalibPars ()

#---------------------

    def loadAllCalibPars (self) :
        """Loads all calibration parameters, if the files are available or set default.
        """
        self.cpars = {}

        for type in self.list_of_clib_types :
            fname = findCalibFile (self.path_to_calib_types, type, self.run) # self.path_to_calib_types + type + '/0-end.data'
            #print 'Load calibpars: ', fname

            cpars_for_type = self.loadCalibParsFromFileOrDefault (fname, type)

            # Special case of 3-d arrays:
            if type == 'center' \
            or type == 'center_corr' \
            or type == 'center_global' : cpars_for_type.shape = (3,4,8)
            #print 'cpars_for_type.shape = ', cpars_for_type.shape
            self.cpars[type] = cpars_for_type

        #=================================
        #self.cpeval = cpe.CSPADCalibParsEvaluated(self)
        #self.cpeval.setCalibParsEvaluated()
        #self.cpeval.printCalibParsEvaluated ('center_global') 
        #=================================

#---------------------

    def loadCalibParsFromFileOrDefault (self, fname, type) :
        """Load parameters of specified type from file or set default.
        """
        
        if fname == None :
            print 'WARNING: CALIBRATION FILE\n', fname, '\nDOES NOT EXIST, WILL USE DEFAULT CONSTANTS.'
            self.cpars_status[type] = 'DEFAULT'
            return self.getCalibParsDefault (type)

        try :
            shape = self.cpars[type].shape # preserve default shape
            self.cpars_status[type] = 'FROM_FILE'
            arr = np.loadtxt (fname)
            arr.shape = shape              # set default shape
            self.cpars[type] = arr
            return self.cpars[type]

        except IOError :
            print 80*'!'
            print 'WARNING: CALIBRATION FILE\n', fname, '\nIS CORRUPTED, WILL USE DEFAULT CONSTANTS.'
            print 80*'!'
            self.cpars_status[type] = 'DEFAULT'
            return self.getCalibParsDefault (type)

#---------------------

    def printCalibFiles (self) :
        """Print the list of calibration files for this object.
        """
        print '\nprintCalibFiles(): List of clib files:'
        msg = ''
        for type in self.list_of_clib_types :
            fname = findCalibFile (self.path_to_calib_types, type, self.run)
            msg += 'Calib type: %s has file: %s\n' % (type.ljust(15), fname)
        print msg

#---------------------

    def printListOfCalibTypes (self) :
        """Print the list of calibration types for this object.
        """
        print 'list_of_clib_types:'
        for type in self.list_of_clib_types : print '   ' + type

#---------------------

    def printCalibPars (self, type=None) :
        """Print all calibration parameters.
        """
        print '\nprintCalibPars(): Calibration parameters:'
        if type==None :
            for type in self.list_of_clib_types :
                print '\nCalibration constants type "' + type + '" with shape' + str(self.cpars[type].shape)
                print self.cpars[type]
        else :
            if type in self.list_of_clib_types :
                print '\nprintCalibParsDefault(): Calibration constants type "' + type + '"' # + '" with shape', self.cpars[type].shape
                print self.cpars[type]
            else :
                msg =  'WARNING: THE REQUESTED TYPE OF CALIBRATION PARS "' + type + \
                       '" IS NOT FOUND IN THE AVAILABLE LIST:\n' + str(self.list_of_clib_types)
                print msg
 
#---------------------

    def printCalibParsStatus (self) :
        """Print status of calibration parameters for all specified files.
        """
        print '\nprintCalibParsStatus(): Status of CSPAD2x2 calibration parameters:'
        for type, status in self.cpars_status.iteritems() :
            print 'Type: %s    Status: %s    Shape: %s' % (type.ljust(16), status.ljust(10), str(self.cpars[type].shape))
 
#---------------------

    def getCalibPars (self, type, run=None) :
        """Returns the numpy array of calibration parameters for specified type and optional run.
        """
        if run!=None : self.run = run

        if not (type in self.list_of_clib_types) :
            msg = 'WARNING: THE REQUESTED TYPE OF CALIBRATION PARS "' + type + \
                   '" IS NOT FOUND IN THE AVAILABLE LIST:\n' + str(self.list_of_clib_types)
            print msg
            return None

        fname = findCalibFile (self.path_to_calib_types, type, self.run)

        #if type == 'center_global' : # for test purpose, in order to call getCalibParsEvaluated (type)
        if type == 'center_global' and fname == None :

            self.cpars_status[type] = 'EVALUATED'

            cpeval = cpe.CSPADCalibParsEvaluated(self)
            #print 'cpeval.getCalibParsEvaluated (type)\n:', cpeval.getCalibParsEvaluated (type)
            return cpeval.getCalibParsEvaluated (type)
            
        return self.loadCalibParsFromFileOrDefault (fname, type)

#---------------------
#---------------------
#---------------------
#---------------------
#---------------------

def findCalibFile (path_to_clib_types, type=None, run=None) :
    """Use the run number, self.path_to_calib_types, and type of the calibration constants.
       From the directory self.path_to_calib_types + '/' + type select the file
       which run range is valid for requested run.
       None is returned if the file is not found.
    """

    err_msg_prefix = 'findCalibFile(): ERROR in findCalibFile(path, type, run): '

    if type==None :
        print  err_msg_prefix + 'type IS NOT SPECIFIED'
        return None

    if run==None :
        print  err_msg_prefix + 'run IS NOT SPECIFIED'
        return None

    if path_to_clib_types[-1] != '/' : path = path_to_clib_types + '/' + type
    else                             : path = path_to_clib_types + type

    run_max = 9999
    calibfname = None

    if not os.path.exists(path) :
        print  'WARNING in findCalibFile(): PATH %s DOES NOT EXIST.' % path
        return calibfname

    flist = os.listdir(path)
    if len(flist) > 1 : flist.sort()

    for fname in flist :

        if fname[-1] == '~' : continue # skip old files with ~(tilde) at the end

        basename = fname.split('.') # Assume: basename[0]='0-end', basename[1]='data' 
        basename_beg, basename_end = basename[0].split('-')

        run_beg = int(basename_beg)
        if basename_end == 'end' :
            run_end = run_max
        else :
            run_end = int(basename_end)

        # Find the last file in the list which run number is in the range
        if run_beg <= run and run <= run_end :
            calibfname = fname

        #print fname, basename[0], run, run_beg, run_end, calibfname

    if calibfname != None : calibfname = path + '/' + calibfname

    #print 'calibfname = ', calibfname
    return calibfname

#----------------------------------------------
# In case someone decides to run this module --
#----------------------------------------------

def main_test() :

    path  = '/reg/d/psdm/xpp/xpptut13/calib/CsPad::CalibV1/XppGon.0:Cspad.0/' 
    run   = 10

    calib = CalibPars(path, run) # , list_of_clib_types=['center', 'tilt', 'pedestals']

    calib.printCalibPars() # prints the default calib pars
    calib.printCalibPars('center_global') 
    print 'Test of getCalibPars("center_global", run):\n', calib.getCalibPars('center_global', run)
    print 'Test of getCalibPars("center", run):\n', calib.getCalibPars('center', run)
    print 'Test of getCalibPars("tilt", run):\n',   calib.getCalibPars('tilt', run)
    calib.printCalibParsStatus()
    calib.printListOfCalibTypes()
    calib.printCalibFiles()

    print 'End of test'

#---------------------

def test_pedestals() :

    path  = '/reg/d/psdm/xpp/xpptut13/calib/CsPad::CalibV1/XppGon.0:Cspad.0/' 
    run   = 10
    calib = CalibPars(path, list_of_clib_types=['center', 'tilt', 'pedestals'])
    arr_pedestals = calib.getCalibPars('pedestals', run)
    print '\narr_pedestals =\n', arr_pedestals
    print '\narr_pedestals.shape =', arr_pedestals.shape

#---------------------

if __name__ == "__main__" :

    if len(sys.argv)==1   :
        main_test()
        print 'For other test(s) use command: python', sys.argv[0], '<test-number=1,2,...>'
    elif sys.argv[1]=='1' : main_test()
    elif sys.argv[1]=='2' : test_pedestals()
    else : print 'Non-expected arguments: sys.argv=', sys.argv

    sys.exit ( 'End of job' )

#----------------------------------------------
