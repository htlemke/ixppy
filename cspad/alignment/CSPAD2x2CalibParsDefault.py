#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module CSPAD2x2CalibParsDefault
#
#------------------------------------------------------------------------

"""
This software was developed for the SIT project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@see RelatedModule

@version $Id: 2013-05-10$

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
import numpy as np

#---------------------
#  Class definition --
#---------------------

class CSPAD2x2CalibParsDefault (object) :
    """Provides access to the CSPAD2x2 calibration parameters through the singleton object
       cspad2x2calibparsdefault.

       This class should not be used by itself; default parameters are different from real.
       It is used in CSPAD2x2CalibPars to get rid of undefined parameters in case of missing calibration files.

       Interface:
       ==========
       Instatiation for singleton cspad2x2calibparsdefault is already done in CSPAD2x2CalibParsDefault, so
       import PyCSPadImage.CSPAD2x2CalibParsDefault as cpd

       Access method:
       center = cpd.cspad2x2calibparsdefault.getCalibParsDefault('center')

       Test methods: 
       cpd.cspad2x2calibparsdefault.printCalibParsDefault()          # For all types
       cpd.cspad2x2calibparsdefault.printCalibParsDefault('center')  # for specified type
       cpd.cspad2x2calibparsdefault.printListOfCalibTypes()
    """

    list_of_clib_types =[
         'center'
        ,'tilt'
        ,'beam_vector'
        ,'common_mode'
        ,'pedestals'
        ,'pixel_status'
        ,'filter'
        ]

#---------------------

    def __init__ (self) :
        """Constructor does not need in input parameters.
           All data implemented calibration types will be initialized.
        """
        self.loadCSPAD2x2CalibParsDefault()

#---------------------

    def loadCSPAD2x2CalibParsDefault (self) :
        """Initialization of all implemented arrays of calibration parameters
        """
        self.defpars = {}

        self.defpars['center'] = np.array(  [[198., 198.],
                                             [ 95., 308.],
                                             [  0.,   0.]])

        self.defpars['tilt']          = np.zeros((2), dtype=np.float32)

        self.defpars['beam_vector']   = np.zeros((3), dtype=np.float32)

        self.defpars['common_mode']   = np.array([1, 100, 30])

        self.defpars['pedestals']     = np.zeros((185, 388, 2), dtype=np.float32)

        self.defpars['pixel_status']  = np.zeros((185, 388, 2), dtype=np.uint16)

        self.defpars['filter']        = np.array([1, 100, 10])

#---------------------

    def printCalibParsDefault (self, partype=None) :
        """Print the calibration prarameters of specified partype or all for dafault.
        """        
        if partype==None :
            for type in self.list_of_clib_types :
                print '\nprintCalibParsDefault(): Calibration constants type "' + type + '"' # + '" with shape', self.cpars[type].shape
                print self.defpars[type]
        else :
            if partype in self.list_of_clib_types :
                print '\nprintCalibParsDefault(): Calibration constants type "' + partype + '"' # + '" with shape', self.cpars[type].shape
                print self.defpars[partype]
            else :
                msg =  'WARNING: THE REQUESTED TYPE OF CALIBRATION PARS "' + partype + \
                       '" IS NOT FOUND IN THE AVAILABLE LIST:\n' + str(self.list_of_clib_types)
                print msg
            
#---------------------

    def printListOfCalibTypes (self) :
        """Print the list of calibration types.
        """
        print '\nprintListOfCalibTypes(): list_of_clib_types:' #, self.list_of_clib_types
        for type in self.list_of_clib_types : print '    ', type

#---------------------

    def getCalibParsDefault (self, type) :
        """Returns calibration parameters of specified type as a numpy array.
        """
        if type in self.list_of_clib_types :
            return self.defpars[type]
        else :
            msg = 'WARNING: THE REQUESTED TYPE OF CALIBRATION PARS "' + type + \
                  '" IS NOT FOUND IN THE AVAILABLE LIST:\n' + str(self.list_of_clib_types)
            print msg
            return None

#---------------------------------------
# Define the default calibration parameters in the singleton object
cspad2x2calibparsdefault = CSPAD2x2CalibParsDefault()

#----------
#-- TEST --
#----------

def main_test() :

    print 'CSPAD2x2CalibParsDefault is enable as a singletone cspad2x2calibparsdefault'    
    cspad2x2calibparsdefault.printCalibParsDefault()
    cspad2x2calibparsdefault.printListOfCalibTypes()
    cspad2x2calibparsdefault.printCalibParsDefault('center')
    print '\nTest of getCalibParsDefault("center"):\n', cspad2x2calibparsdefault.getCalibParsDefault('center')

if __name__ == "__main__" :
    main_test()
    sys.exit ( 'End of job' )

#----------------------------------------------
