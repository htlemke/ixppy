#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module CSPAD2x2CalibPars...
#
#------------------------------------------------------------------------

"""
This software was developed for the SIT project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@see RelatedModule

@version $Id: 2013-05-10$

@author Mikhail S. Dubrovin
"""

#---------------------
import sys
import os

import numpy as np
import CSPAD2x2CalibParsDefault as cpd
from   CalibPars import findCalibFile

#---------------------
#  Class definition --
#---------------------

class CSPAD2x2CalibPars (object) :
    """This class provides access to the CSPAD2x2 calibration parameters.

       Interface
       =========
       Regular instantiation:
           ALL parameters are OPTIONAL NAMED parameters;
           path  = '/reg/d/psdm/mec/mec73313/calib/CsPad2x2::CalibV1/MecTargetChamber.0:Cspad2x2.1/'
           run   = 123
           list_of_clib_types = ['center', 'tilt', 'pedestals'] # optional list of used types
           calib = CSPAD2x2CalibPars(path, run, list_of_clib_types)

       Other option for instantiation:
           calib    = CSPAD2x2CalibPars() # list_of_clib_types - is an optional
           run      = 123  - is an optional, named
           calibdir = '/reg/d/psdm/mec/mec73313/calib'
           group    = 'CsPad2x2::CalibV1'
           source   = 'MecTargetChamber.0:Cspad2x2.1'
           calib.setCalibPars (run, calibdir, group, source)

       Get array of calibration parameters for specified type and run number:
           type = 'center'
           arr  = calib.getCalibPars (type[,run])
    """
    #enum_status = {0:'DEFAULT', 1:'FROM_FILE'}

    list_of_clib_types_total = cpd.cspad2x2calibparsdefault.list_of_clib_types

#---------------------

    def __init__ (self, path=None, run=None, list_of_clib_types=None) :
        """Class constructor:
           The path and run need to be defined to instantiate the object and load correct set of parameters.
           If path or run is omitted, default parameters will be used. 
           Run number can be omitted here and passed later in getCalibPars(type, run)
           list_of_clib_types - optional parameter for optimization of time; only types from the list are loaded.
           If the list_of_clib_types is omitted, all parameters will be loaded and used.
        """    
        self.cpars        = {}
        self.cpars_status = {} # 'DEFAULT', 'FROM_FILE'

        self.run = run
        self.path_to_calib_types = path

        self.setListOfCalibTypes(list_of_clib_types)
        
        self.setCalibParsDefault()

        #if path!=None and run!=None :
        #    self.setCalibParsForPath(path, run)

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
        return cpd.cspad2x2calibparsdefault.getCalibParsDefault (type)

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
                      run      = None,   # 123
                      calibdir = None,   # '/reg/d/psdm/mec/mec73313/calib'
                      group    = None,   # 'CsPad2x2::CalibV1'
                      source   = None) : # 'MecTargetChamber.0:Cspad2x2.1'
        """Set path and run for calibration parameters.
        """
        
        self.setCalibParsForPath (calibdir + '/' + group + '/' + source, run)

#---------------------

    def setCalibParsForPath (self,
                             path, # '/reg/d/psdm/mec/mec73313/calib/CsPad2x2::CalibV1/MecTargetChamber.0:Cspad2x2.1'
                             run  = None, # 123
                             ) :
        """ Set path and run for calibration parameters.
        """
        self.path_to_calib_types = path
        if run!=None  : self.run = run

        # self.loadAllCalibPars ()

#---------------------

    def loadAllCalibPars (self) :
        """Loads all calibration parameters, if the files are available or set default.
        """
        self.cpars = {}

        for type in self.list_of_clib_types :
            fname = findCalibFile (self.path_to_calib_types, type, self.run)
            #print 'Load calibpars: ', fname

            cpars_for_type = self.loadCalibParsFromFileOrDefault (fname, type)

            print 'cpars_for_type.shape = ', cpars_for_type.shape
        
            # Special case of array shapes:
            if type == 'pedestals' \
            or type == 'pixel_status' : cpars_for_type.shape = (185,388,2)

            self.cpars[type] = cpars_for_type

#---------------------

    def loadCalibParsFromFileOrDefault (self, fname, type) :
        """Load parameters of specified type from file or set default.
        """

        if fname == None :
            self.cpars_status[type] = 'DEFAULT'
            #print 'WARNING: CALIBRATION FILE ', fname, '\nDOES NOT EXIST, WILL USE DEFAULT CONSTANTS.'
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
        print '\nprintListOfCalibTypes(): List of clib types:'
        for type in self.list_of_clib_types : print '   ' + type

#---------------------

    def printCalibPars (self, type=None) :
        """Print all calibration parameters.
        """
        print '\nprintCalibPars(): Calibration parameters:'
        if type==None :
            for t in self.list_of_clib_types :
                print '\nCalibration constants type "' + t + '" with shape' + str(self.cpars[t].shape)
                print self.cpars[t]

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
            print 'Type: %s    Status: %s    Shape: %s' % (type.ljust(12), status.ljust(10), str(self.cpars[type].shape))

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
        return self.loadCalibParsFromFileOrDefault (fname, type)

#---------------------
#---------------------
#---------------------
#------- TEST --------
#---------------------
#---------------------
#---------------------

def main_test() :

    #path = '/reg/d/psdm/mec/mec73313/calib/CsPad2x2::CalibV1/MecTargetChamber.0:Cspad2x2.1/'
    #path = '/reg/neh/home1/dubrovin/LCLS/CSPad2x2Alignment/calib-cspad2x2-02-2013-02-13/'
    path  = '/reg/d/psdm/xpp/xpptut13/calib/CsPad2x2::CalibV1/XppGon.0:Cspad2x2.0/'
    run   = 180

    #calib = CSPAD2x2CalibPars() # Sets all default calibration parameters
    #calib = CSPAD2x2CalibPars(list_of_clib_types=['center', 'tilt', 'pedestals']) # Sets default calibration parameters
    #calib.setCalibParsForPath(run, path)                           
    #calib = CSPAD2x2CalibPars(path, run)
    calib = CSPAD2x2CalibPars(path, run) #, ['center', 'tilt', 'pedestals'])

    calib.printCalibPars()
    calib.printCalibPars ('center')
    print 'Test of getCalibPars("center", run):\n', calib.getCalibPars('center', run)
    print 'Test of getCalibPars("tilt", run):\n',   calib.getCalibPars('tilt', run)
    calib.printCalibParsStatus()
    calib.printListOfCalibTypes()
    calib.printCalibFiles()


def test_pedestals() :

    path  = '/reg/d/psdm/xpp/xpptut13/calib/CsPad2x2::CalibV1/XppGon.0:Cspad2x2.1/'
    calib = CSPAD2x2CalibPars(path, list_of_clib_types=['center', 'tilt', 'pedestals'])
    run   = 180
    arr_pedestals = calib.getCalibPars('pedestals', run)
    print '\narr_pedestals =\n', arr_pedestals
    print '\narr_pedestals.shape =', arr_pedestals.shape


if __name__ == "__main__" :

    if len(sys.argv)==1   :
        main_test()
        print 'For other test(s) use command: python', sys.argv[0], '<test-number=1,2,...>'
    elif sys.argv[1]=='1' : main_test()
    elif sys.argv[1]=='2' : test_pedestals()
    else : print 'Non-expected arguments: sys.argv=', sys.argv

    sys.exit ( 'End of job' )

#----------------------------------------------
