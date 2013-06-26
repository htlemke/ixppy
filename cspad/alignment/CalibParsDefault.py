#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module CalibParsDefault...
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

#---------------------
#  Class definition --
#---------------------

class CalibParsDefault (object) :
    """This class provides access to the calibration parameters
    """

#---------------------

    def __init__ (self) :
        """Constructor"""

        self.list_of_clib_types =[  'center'
                                   ,'center_corr'  
                                   ,'marg_gap_shift' 
                                   ,'offset'
                                   ,'offset_corr'
                                   ,'rotation'
                                   ,'tilt'
                                   ,'quad_rotation'
                                   ,'quad_tilt'
                                   ,'common_mode'
                                   ,'pedestals'
                                   ,'filter'
                                   ,'pixel_status'
                                   ,'center_global'  
                                   ,'rotation_index'
                                   ]

        self.loadCalibParsDefault()

#---------------------

    def loadCalibParsDefault (self) :

        self.defpars = {}

        self.defpars['center'] = np.array(
                    [[[198.,  198.,  310.,   98.,  627.,  628.,  711.,  498.],
                      [198.,  198.,  310.,   98.,  627.,  628.,  711.,  498.],
                      [198.,  198.,  310.,   98.,  627.,  628.,  711.,  498.],
                      [198.,  198.,  310.,   98.,  627.,  628.,  711.,  498.]],
        
                     [[307.,   95.,  625.,  625.,  515.,  727.,  198.,  199.],
                      [307.,   95.,  625.,  625.,  515.,  727.,  198.,  199.],
                      [307.,   95.,  625.,  625.,  515.,  727.,  198.,  199.],
                      [307.,   95.,  625.,  625.,  515.,  727.,  198.,  199.]],
       
                     [[  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
                      [  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
                      [  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
                      [  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.]]])

        self.defpars['center_corr'] = np.array(
                    [[[  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
                      [  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
                      [  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
                      [  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.]],
                     [[  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
                      [  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
                      [  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
                      [  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.]],
                     [[  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
                      [  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
                      [  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.],
                      [  0.,    0.,    0.,    0.,    0.,    0.,    0.,    0.]]])

        self.defpars['marg_gap_shift'] = np.array(
                    [[ 15.,  40.,   0.,  38.],
                     [ 15.,  40.,   0.,  38.],
                     [  0.,   0.,   0.,   0.]])

        self.defpars['offset'] = np.array(
                    [[   0.,    0.,  834.,  834.],
                     [   0.,  834.,  834.,    0.],
                     [   0.,    0.,    0.,    0.]])

        self.defpars['offset_corr'] = np.array(
                    [[   0.,    0.,    0.,    0.],
                     [   0.,    0.,    0.,    0.],
                     [   0.,    0.,    0.,    0.]])

        self.defpars['rotation'] = np.array(
                    [[   0.,    0.,  270.,  270.,  180.,  180.,  270.,  270.],
                     [   0.,    0.,  270.,  270.,  180.,  180.,  270.,  270.],
                     [   0.,    0.,  270.,  270.,  180.,  180.,  270.,  270.],
                     [   0.,    0.,  270.,  270.,  180.,  180.,  270.,  270.]])

        self.defpars['tilt'] = np.array(
                    [[0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],  
                     [0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],  
                     [0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.],  
                     [0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.]])

        self.defpars['quad_rotation'] = np.array([ 180.,   90.,    0.,  270.])

        self.defpars['quad_tilt']     = np.array([   0.,    0.,    0.,    0.])

        self.defpars['common_mode']   = np.array([   1, 100, 30])

        self.defpars['filter']        = np.array([   1, 100, 10])

        self.defpars['pedestals']     = np.zeros((5920, 388), dtype=np.float32) # SHAPE: (5920, 388)

        self.defpars['pixel_status']  = np.zeros((5920, 388), dtype=np.uint16) # SHAPE: (5920, 388)

        self.defpars['center_global'] = np.array(
           [[[ 473.38,  685.26,  155.01,  154.08,  266.81,   53.95,  583.04,  582.15],  
             [ 989.30,  987.12, 1096.93,  884.11, 1413.16, 1414.94, 1500.83, 1288.02],  
             [1142.59,  930.23, 1459.44, 1460.67, 1347.57, 1559.93, 1032.27, 1033.44],  
             [ 626.78,  627.42,  516.03,  729.15,  198.28,  198.01,  115.31,  327.66]],  

            [[1028.07, 1026.28, 1139.46,  926.91, 1456.78, 1457.35, 1539.71, 1327.89],  
             [1180.51,  967.36, 1497.74, 1498.54, 1385.08, 1598.19, 1069.65, 1069.93],  
             [ 664.89,  666.83,  553.60,  765.91,  237.53,  236.06,  152.17,  365.47],  
             [ 510.38,  722.95,  193.33,  193.41,  308.04,   95.25,  625.28,  624.14]],  
       
            [[     0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.],
             [     0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.],
             [     0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.],
             [     0.,      0.,      0.,      0.,      0.,      0.,      0.,      0.]]])

        self.defpars['rotation_index'] = np.array(
                 [[   0,   0,   3,   3,   2,   2,   3,   3],
                  [   3,   3,   2,   2,   1,   1,   2,   2],
                  [   2,   2,   1,   1,   0,   0,   1,   1],
                  [   1,   1,   0,   0,   3,   3,   0,   0]])

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
                print  'WARNING: THE REQUESTED TYPE OF CALIBRATION PARS "', partype, \
                       '" IS NOT FOUND IN THE AVAILABLE LIST:\n', self.list_of_clib_types
            
#---------------------

    def printListOfCalibTypes (self) :
        print '\nprintListOfCalibTypes(): list_of_clib_types:' #, self.list_of_clib_types
        for type in self.list_of_clib_types : print '    ', type

#---------------------

    def getCalibParsDefault (self, type) :

        if type in self.list_of_clib_types :
            return self.defpars[type]
        else :
            print  'WARNING: THE REQUESTED TYPE OF CALIBRATION PARS "', type, \
                   '" IS NOT FOUND IN THE AVAILABLE LIST:\n', self.list_of_clib_types
            return None

#---------------------------------------

calibparsdefault = CalibParsDefault()

#----------------------------------------------
# In case someone decides to run this module --
#----------------------------------------------

def main() :

    calibparsdefault.printCalibParsDefault()
    calibparsdefault.printListOfCalibTypes()
    print 'End of test'

if __name__ == "__main__" :

    main()
    sys.exit ( 'End of job' )

#----------------------------------------------
