#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module GlobalMethods...
#
#------------------------------------------------------------------------

""" Define a bunch of global methods

This software was developed for the SIT project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@see RelatedModule

@version $Id: template!python!py 4 2008-10-08 19:27:36Z salnikov $

@author Mikhail S. Dubrovin
"""


#------------------------------
#  Module's version from CVS --
#------------------------------
__version__ = "$Revision: 4 $"
# $Source$

#--------------------------------
#  Imports of standard modules --
#--------------------------------
#import sys
import os
import sys
import string as st
import time
import numpy as np
#-----------------------------
# Imports for other modules --
#-----------------------------
#from PkgPackage.PkgModule import PkgClass

#import ConfigParameters as cp

#------------------------
# Exported definitions --
#------------------------

#----------------------------------

def get_item_last_name(dsname):
    """Returns the last part of the full item name (after last slash)"""
    path,name = os.path.split(str(dsname))
    return name

#----------------------------------

def get_item_path_to_last_name(dsname):
    """Returns the path to the last part of the item name"""
    path,name = os.path.split(str(dsname))
    return path

#----------------------------------

def get_item_path_and_last_name(dsname):
    """Returns the path and last part of the full item name"""
    path,name = os.path.split(str(dsname))
    return path, name

#----------------------------------

def get_item_second_to_last_name(dsname):
    """Returns the 2nd to last part of the full item name"""
    path1,name1 = os.path.split(str(dsname))
    path2,name2 = os.path.split(str(path1))
    return name2 

#----------------------------------

def get_item_third_to_last_name(dsname):
    """Returns the 3nd to last part of the full item name"""
    path1,name1 = os.path.split(str(dsname))
    path2,name2 = os.path.split(str(path1))
    path3,name3 = os.path.split(str(path2))
    str(name3)
    return name3 

#----------------------------------

def get_item_last_three_names(dsname):
    """Returns the 3nd to last part of the full item name"""
    path1,name1 = os.path.split(str(dsname))
    path2,name2 = os.path.split(str(path1))
    path3,name3 = os.path.split(str(path2))
    #str(name3)
    return name1, name2, name3 

#----------------------------------

def get_item_name_for_title(dsname):
    """Returns the last 3 parts of the full item name (after last slashes)"""
    path1,name1 = os.path.split(str(dsname))
    path2,name2 = os.path.split(str(path1))
    path3,name3 = os.path.split(str(path2))

    return name3 + '/' + name2 + '/' + name1

#----------------------------------

def ImageIsInTheName(dsname):
    path1,name1 = os.path.split(str(dsname))
    path2,name2 = os.path.split(str(path1))
    path3,name3 = os.path.split(str(path2))

    imageIsInTheName = False
    if   name1 == 'image' : imageIsInTheName = True
    elif name1 == 'data' and name3[0:16] == 'Princeton::Frame' : imageIsInTheName = True

    print 'imageIsInTheName :',
    print '       last name:', name1
    print '2nd to last name:', name2
    print '3rd to last name:', name3
    print 'name3[0:16]', name3[0:16]
    print 'imageIsInTheName returns:', imageIsInTheName

    return imageIsInTheName

#----------------------------------

def CSpadIsInTheName(dsname):
    path1,name1 = os.path.split(str(dsname))
    path2,name2 = os.path.split(str(path1))
    path3,name3 = os.path.split(str(path2))

    #print '       last name:', name1
    #print '2nd to last name:', name2
    #print '3rd to last name:', name3
    #print 'name3[0:5]', name3[0:5]

    cspadIsInTheName = False
    if name3[0:5] == 'CsPad' and name1 == 'data' : cspadIsInTheName = True
    #print 'cspadIsInTheName :', cspadIsInTheName

    return cspadIsInTheName

#----------------------------------

def CSpadElementIsInTheName(dsname):
    path1,name1 = os.path.split(str(dsname))
    path2,name2 = os.path.split(str(path1))
    path3,name3 = os.path.split(str(path2))
    if name3[0:14] == 'CsPad::Element' and name1 == 'data' : return True
    else                                                   : return False

#----------------------------------

def CSpad2x2ElementIsInTheName(dsname):
    path1,name1 = os.path.split(str(dsname))
    path2,name2 = os.path.split(str(path1))
    path3,name3 = os.path.split(str(path2))
    if   name3[0:8]  == 'CsPad2x2' and name1 == 'data'           : return True
    elif name3[0:18] == 'CsPad::MiniElement' and name1 == 'data' : return True
    else                                                         : return False



#----------------------------------

#def CSpadDatasetIsChecked():
#    for dsname in cp.confpars.list_of_checked_item_names :
#        if CSpadIsInTheName(dsname) : return True
#    return False

#----------------------------------

#def ImageDatasetIsChecked():
#    for dsname in cp.confpars.list_of_checked_item_names :
#        if ImageIsInTheName(dsname) : return True
#    return False

#----------------------------------

#def WaveformDatasetIsChecked():
#    for dsname in cp.confpars.list_of_checked_item_names :
#        if get_item_last_name(dsname) == 'waveforms' : return True
#    return False

#----------------------------------

#def CorrelationDatasetIsChecked():
#    for dsname in cp.confpars.list_of_checked_item_names :
#        item_last_name = get_item_last_name(dsname)
#        if item_last_name == 'time': return True
#        if item_last_name == 'data' and (not CSpadIsInTheName(dsname)): return True
#    return False

#----------------------------------

def getPatternEndsInTheString(symbolic_string, pattern='CalibCycle:'):

    start = st.find(symbolic_string, pattern)
    pattern_length = len(pattern)
    if start<0 : 
        return -1, -1, False
    else :
        return start, start+pattern_length, True
    
#----------------------------------

def CalibCycleIsInThePath(path_and_name):

    pattern = 'CalibCycle:'
    start = st.find(path_and_name,pattern)
    if start>0 : 
        #name = path_and_name[start:start+11] 
        #numb = path_and_name[start+11:start+11+4]
        #print 'name,numb=',name,numb
        return True
    else :
        return False
    
#----------------------------------

#def CalibCycleDatasetIsChecked():
#    for dsname in cp.confpars.list_of_checked_item_names :
#        if CalibCycleIsInThePath(dsname):
#            print 'CalibCycleDatasetIsChecked(): True'
#            return True
#
#    print 'CalibCycleDatasetIsChecked(): False'
#    return False

#----------------------------------

def getListOfFilesInDir(dirname) :
    return os.listdir(dirname)

#----------------------------------

def printListOfFilesInDir(dirname) :
    dirList = getListOfFilesInDir(dirname)
    print 'List of files in the dir.', dirname
    for name in dirList :
        print name,
    print '\n'

#----------------------------------

def getListOfFiles(dirname, extension='h5') :
    print """Returns the list of files in the specified directory with requested extension"""
    
    if not os.path.exists(dirname) :
        print  'WARNING: THE SPECIFIED DIRECTORY "',dirname,'" DOES NOT EXIST.'
        return None

    dot_extension = '.' + extension
    list_of_files = []
    for fname in os.listdir(dirname) :
        name, ext = os.path.splitext(fname)
        #print fname, ext
        if ext == dot_extension :
            list_of_files.append( dirname + '/' + fname )

    if len(list_of_files) == 0 : return None
    else :                       return list_of_files

#----------------------------------

def saveNumpyArrayInFile(arr, fname='nparray.txt', format='%f') : # format='%f'
    print """Save numpy array in file """, fname
    np.savetxt(fname, arr, fmt=format)

#----------------------------------

def getNumpyArrayFromFile(fname='nparray.txt', datatype=np.float32) : # np.int16, np.float16, np.float32
    print """Load numpy array from file """, fname
    return np.loadtxt(fname, dtype=datatype)

#----------------------------------

def getCSPadArrayFromFile(fname, dtype=np.float32, shape = (32, 185, 388)) : # raw shape is (5920, 388)
    try :
        arr = np.loadtxt(fname, dtype)
        print 'Load array from file:', fname
        #print 'Input  arr.shape=', arr.shape
        #print 'Output arr.shape=', shape
        arr.shape = shape
        return arr

    except IOError:
        print 'IOError: CAN NOT OPEN FILE:', fname
        sys.exit ( "THIS PROBLEM NEEDS TO BE FIXED BEFORE YOU CAN CONTINUE..." )
        #return None

#----------------------------------
# Test

def main() :    
    for fname in getListOfFiles('/reg/d/psdm/cxi/cxii0211/hdf5', 'h5') : print fname

#----------------------------------

if __name__ == "__main__" :

    main()
    sys.exit ( "Module is not supposed to be run as main module" )

#----------------------------------
