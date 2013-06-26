#--------------------------------------------------------------------------
# File and Version Information:
#  $Id$
#
# Description:
#  Module GlobalGraphics...
#
#------------------------------------------------------------------------

"""This module provides access to the calibration parameters

This software was developed for the SIT project.  If you use all or 
part of it, please give an appropriate acknowledgment.

@see RelatedModule

@version $Id: 2008-09-28$

@author Mikhail S. Dubrovin
"""

#------------------------------
#  Module's version from CVS --
#------------------------------
__version__ = "$Revision: 4 $"
# $Source$

#--------------------------------

import sys
import numpy as np

import matplotlib.pyplot as plt

#---------------------
#class CSPadImageProducer (object) :
#    """This is an empty class"""
#---------------------
#    def __init__ (self) :
#        print 'CSPadImageProducer __init__'
#--------------------------------

def getArrangedImage() :
    arr = np.arange(2400)
    arr.shape = (40,60)
    return arr

#--------------------------------
def getRandomImage() :
    mu, sigma = 200, 25
    arr = mu + sigma*np.random.standard_normal(size=2400)
    arr.shape = (40,60)
    return arr

#--------------------------------

def plotHistogram(arr,range=(0,500),figsize=(5,5)) :
    fig  = plt.figure(figsize=figsize, dpi=80, facecolor='w',edgecolor='w',frameon=True)
    plt.hist(arr.flatten(), bins=100, range=range)

#--------------------------------

def plotSpectrum(arr,range=(0,500),figsize=(5,5)) :
    plotHistogram(arr,range,figsize)

#--------------------------------

def plotImage(arr,range=(0,500),figsize=(12,5)) :
    fig  = plt.figure(figsize=figsize, dpi=80, facecolor='w',edgecolor='w',frameon=True)
    axes1 = plt.imshow(arr, interpolation='nearest', aspect='auto') #,extent=self.XYRange, origin='lower'
    colb1 = fig.colorbar(axes1, pad=0.01, fraction=0.1, shrink=1, aspect=20)
    axes1.set_clim(range[0],range[1])

#--------------------------------

def plotImageLarge(arr,range=(0,500),figsize=(12,10)) :
    fig  = plt.figure(figsize=figsize, dpi=80, facecolor='w',edgecolor='w',frameon=True)
    axim = fig.add_axes([0.04,  0.02, 0.99, 0.96])
    axes1 = axim.imshow(arr, interpolation='nearest', aspect='auto') #,extent=self.XYRange, origin='lower'
    colb1 = fig.colorbar(axes1, pad=0.005, fraction=0.09, shrink=1, aspect=40)
    axes1.set_clim(range[0],range[1])

#--------------------------------

def plotImageAndSpectrum(arr,range=(0,500)) :
    fig  = plt.figure(figsize=(15,5), dpi=80, facecolor='w',edgecolor='w',frameon=True)
    fig.canvas.set_window_title('Image And Spectrum ' + u'\u03C6')

    ax1   = plt.subplot2grid((10,10), (0,4), rowspan=10, colspan=6)
    axim1 = ax1.imshow(arr, interpolation='nearest', aspect='auto') # , origin='lower' 
    colb1 = fig.colorbar(axim1, pad=0.01, fraction=0.1, shrink=1.00, aspect=20)
    axim1.set_clim(range[0],range[1])
    plt.title('Image',color='b',fontsize=20)

    ax2   = plt.subplot2grid((10,10), (0,0), rowspan=10, colspan=4)
    ax2.hist(arr.flatten(), bins=100, range=range)
    plt.title('Spectrum',color='r',fontsize=20)
    plt.xlabel('Bins')
    plt.ylabel('Stat') #u'\u03C6'

#--------------------------------

def move(x0=200,y0=100) :
    #plt.get_current_fig_manager().window.move(x0, y0)
    move_str = '+' + str(x0) + '+' + str(y0)
    plt.get_current_fig_manager().window.geometry(move_str)
    #plt.get_current_fig_manager().window.geometry("+50+50")
    pass

#--------------------------------

def show() :
    plt.show()
    #file.close()

#----------------------------------------------

def main() :

    arr = getRandomImage()
    #plotImage(arr,range=(100,300))
    plotImageAndSpectrum(arr,range=(100,300))
    #plotHistogram(arr,range=(100,300),figsize=(10,5))
    #plotImage(arr,range=(100,300),figsize=(10,10))
    move(200,100)
    show()

#--------------------------------

if __name__ == "__main__" :

    main()
    sys.exit ( 'End of test.' )

#----------------------------------------------
