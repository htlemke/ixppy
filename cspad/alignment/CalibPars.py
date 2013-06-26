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
import CalibParsDefault   as cpdef
import CalibParsEvaluated as cpe

#---------------------
#  Class definition --
#---------------------

class CalibPars (object) :
    """This class provides access to the calibration parameters
    """

#---------------------

    def __init__ (self) :

        self.run = None # 10
        self.path_to_calib_types = None #'/reg/d/psdm/CXI/cxi35711/calib'+'/'+'CsPad::CalibV1'+'/'+'CxiDs1.0:Cspad.0'
        
        self.list_of_clib_types =[  'center'
                                   ,'center_corr'  
                                   ,'marg_gap_shift' 
                                   ,'offset'
                                   ,'offset_corr'
                                   ,'rotation'
                                   ,'tilt'
                                   ,'quad_rotation'
                                   ,'quad_tilt'
                                   #,'common_mode'
                                   #,'filter'
                                   #,'pedestals'
                                   #,'pixel_status'
                                   ]

        self.setCalibParsDefault()
        #self.setCalibPars(10, '/reg/d/psdm/CXI/cxi35711/calib', 'CsPad::CalibV1', 'CxiDs1.0:Cspad.0')

#---------------------

    def setCalibParsDefault (self) :

        self.cpars = {}
        #print 'Set default calibration parameters'
        for type in self.list_of_clib_types :
            self.cpars[type] = self.getCalibParsDefault (type)

#---------------------

    def getCalibParsDefault (self, type) :
        return cpdef.calibparsdefault.getCalibParsDefault (type)

#---------------------

    def setRun ( self, run = 0 ) :
        self.run = run

        #print 'Load the calibration parameters for run ', self.run
        self.setCalibPars()

#---------------------

    def setCalibPars (self,
                      run      = None, # 1
                      calibdir = '/reg/d/psdm/CXI/cxi35711/calib',
                      group    = 'CsPad::CalibV1',
                      source   = 'CxiDs1.0:Cspad.0') :
        """ Set calibration parameters for specified input pars"""

        self.run = run
        self.path_to_calib_types = calibdir + '/' + group + '/' + source + '/'
        self.loadAllCalibPars ()

#---------------------

    def setCalibParsForPath (self,
                      run      = None, # 1
                      path     = '/reg/d/psdm/CXI/cxi35711/calib/CsPad::CalibV1/CxiDs1.0:Cspad.0' ) :
        """ Set calibration parameters for specified input pars"""

        self.run = run
        self.path_to_calib_types = path + '/'
        self.loadAllCalibPars ()

#---------------------

    def loadAllCalibPars (self) :

        self.cpars = {}

        #print 'Load the calibration parameters for run ', self.run

        for type in self.list_of_clib_types :
            fname = self.findCalibFile (self.run, type) # self.path_to_calib_types + type + '/0-end.data'
            #print 'Load calibpars: ', fname

            calibpars = self.loadCalibParsFromFileOrDefault (fname, type)
            if type == 'center' or type == 'center_corr' : calibpars.shape = (3,4,8)
            #print 'calibpars.shape = ', calibpars.shape
            self.cpars[type] = calibpars

        #=================================
        cpe.cpeval.setCalibParsEvaluated()
        #cpe.cpeval.printCalibParsEvaluated ('center_global') 
        #=================================

#---------------------

    def loadCalibParsFromFileOrDefault (self, fname, type) :

        if fname == None :
            return self.getCalibParsDefault (type)

        try :
            return np.loadtxt (fname)

        except IOError :
            print 80*'!'
            print 'WARNING: CALIBRATION FILE\n', fname, '\nDOES NOT EXIST OR CORRUPTED, WILL USE DEFAULT CONSTANTS.'
            print 80*'!'
            return self.getCalibParsDefault (type)

#---------------------

    def printCalibFiles (self) :

        for type in self.list_of_clib_types :
            fname = self.findCalibFile (self.run, type) # self.path_to_calib_types + type + '/0-end.data'
            print 'Print calibpars fname: ', fname

#---------------------

    def printListOfCalibTypes (self) :
        print 'list_of_clib_types:', self.list_of_clib_types

#---------------------

    def printCalibPars (self) :

        for type in self.list_of_clib_types :
            print '\nCalibration constants type "' + type + '" with shape', self.cpars[type].shape
            print self.cpars[type]
            
#---------------------

    def getCalibPars (self, type) :

        if type in self.list_of_clib_types :
            return self.cpars[type]
        else :
            print  'WARNING: THE REQUESTED TYPE OF CALIBRATION PARS "', type, \
                   '" IS NOT FOUND IN THE AVAILABLE LIST:\n', self.list_of_clib_types
            return None

#---------------------

    def findCalibFile (self, run=0, type='offset_corr') :
        """Use the run number, self.path_to_calib_types, and type of the calibration constants.
           From the directory self.path_to_calib_types + '/' + type select the file
           which run range is valid for requested run.
           None is returned if the file is not found.
        """

        path = self.path_to_calib_types + type

        self.run_max = 9999
        self.calibfname = None

        if not os.path.exists(path) :
            print  'WARNING: THE SPECIFIED DIRECTORY "',path,'" DOES NOT EXIST.'
            return self.calibfname

        flist = os.listdir(path)
        if len(flist) > 1 : flist.sort()

        for fname in flist :
            basename = fname.split('.') # Assume: basename[0]='0-end', basename[1]='data' 
            basename_beg, basename_end = basename[0].split('-')

            self.run_beg = int(basename_beg)
            if basename_end == 'end' :
                self.run_end = self.run_max
            else :
                self.run_end = int(basename_end)

            # Find the last file in the list which run number is in the range
            if self.run_beg <= run and run <= self.run_end :
                self.calibfname = fname

            #print fname, basename[0], run, self.run_beg, self.run_end, self.calibfname

        if self.calibfname != None : self.calibfname = path + '/' + self.calibfname

        #print 'self.calibfname = ', self.calibfname
        return self.calibfname

#---------------------------------------

calibpars = CalibPars() # Sets default calibration parameters.

#----------------------------------------------
# In case someone decides to run this module --
#----------------------------------------------

def main() :

    calibpars.printListOfCalibTypes()
    calibpars.printCalibPars() # prints the default calib pars

    #calibpars.setCalibPars(10, '/reg/d/psdm/CXI/cxi35711/calib', 'CsPad::CalibV1', 'CxiDs1.0:Cspad.0')
    calibpars.setCalibParsForPath (run=10, path='/reg/d/psdm/CXI/cxi35711/calib/CsPad::CalibV1/CxiDs1.0:Cspad.0')
    calibpars.printCalibPars()
    cpe.cpeval.printCalibParsEvaluatedAll() 

    #calibpars.printCalibFiles()
    #calibpars.findCalibFile(999)
    #print calibpars.getCalibPars('offset')
    #print calibpars.getCalibPars('XxX')
    #cpdef.calibparsdefault.printCalibParsDefault()

    print 'End of test'

if __name__ == "__main__" :

    main()
    sys.exit ( 'End of job' )

#----------------------------------------------
