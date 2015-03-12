import numpy as np
import scipy as sp
import scipy.signal as signal
import matplotlib
import pylab as pl
import time
import types
import numpy.ma as ma
import os
from toolsVecAndMat import *
from toolsDistrAndHist import *

from matplotlib.lines import Line2D
from matplotlib.widgets import Widget
from copy import copy 
from matplotlib import pyplot as plt
from distutils.version import LooseVersion
from matplotlib.ticker import MaxNLocator


if LooseVersion(matplotlib.__version__)<'1.2':
  from matplotlib.nxutils import points_inside_poly
else:
  from matplotlib.path import Path as mPath
  

def nfigure(name="noname",figsize=None,**figprops):
	try:
		fig_names = [x.canvas.manager.window.get_title()
				for x in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
	except:
          try:
		fig_names = [x.canvas.manager.window.wm_title()
				for x in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]
          except:
		fig_names = [x.canvas.manager.window.windowTitle()
				for x in matplotlib._pylab_helpers.Gcf.get_all_fig_managers()]

	#import code; code.interact(local=locals())

	n=0
	found=0
	for tname in fig_names:
		n+=1
		if tname == name:
			fig=matplotlib._pylab_helpers.Gcf.get_all_fig_managers()[n-1]
			matplotlib._pylab_helpers.Gcf.set_active(fig)
			fig = pl.gcf()
			found = 1

	if not found==1:
		print 'Created new figure %s'  % (name)
		if (figsize is None):
			fig = pl.figure(**figprops)
		else:
			fig = pl.figure(figsize=figsize,**figprops)
		if name is not None:
		  fig.canvas.set_window_title(name)
#	figure(fig.number)
	return fig

def draw_verticalline(pos=0,linespec='k'):
  ah = pl.gca()
  yl = ah.get_ylim()
  pl.plot(pos*np.ones(2),yl,linespec)

def draw_horizontalline(pos=0,linespec='k'):
  ah = pl.gca()
  xl = ah.get_xlim()
  pl.plot(xl,pos*np.ones(2),linespec)

def mfreqz(b,a=1):
    w,h = signal.freqz(b,a)
    h_dB = 20 * np.log10 (abs(h))
    pl.subplot(211)
    pl.plot(w/max(w),h_dB)
    pl.ylim(-150, 5)
    pl.ylabel('Magnitude (db)') 
    pl.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    pl.title(r'Frequency response')
    pl.subplot(212)
    h_Phase = np.unwrap(np.arctan2(np.imag(h),np.real(h)))
    pl.plot(w/max(w),h_Phase)
    pl.ylabel('Phase (radians)') 
    pl.xlabel(r'Normalized Frequency (x$\pi$rad/sample)')
    pl.title(r'Phase response')
    pl.subplots_adjust(hspace=0.5)

def impz(b,a=1):
        impulse = np.repeat(0.,50); impulse[0] =1.
        x = np.arange(0,50)
        response = signal.lfilter(b,a,impulse)
        pl.subplot(211)
        pl.stem(x, response)
        pl.ylabel('Amplitude') 
        pl.xlabel(r'n (samples)')
        pl.title(r'Impulse response')
        pl.subplot(212)
        step = np.cumsum(response)
        pl.stem(x, step)
        pl.ylabel('Amplitude') 
        pl.xlabel(r'n (samples)')
        pl.title(r'Step response')
        pl.subplots_adjust(hspace=0.5)


def clim_std(stdfac=1,axh=None):
  if not axh:
    axh = pl.gca()
  ihandles = axh.findobj(matplotlib.image.AxesImage)
  ih = ihandles[-1]
  im = ih.get_array()
  imean = np.median(im.ravel())
  istd = mad(im.ravel())
  ih.set_clim([imean-stdfac*istd,imean+stdfac*istd])

#TODO section

def getSpanCoordinates(direction='horizontal',axh=None,fig=None,data=None):
  """Tool for selecting a span, functionality similar to ginput. Finish with right mouse button."""
  if not axh:
    axh = pl.gca()
  if not fig: fig = pl.gcf()
  class ROI:
    def __init__(self,fig,axh,direction):
      self.fig = fig
      self.axh = axh
      self.lims = []
      self.boxh = []
      self.finished = False
      self.direction = direction

    def coo(self,tmin,tmax):
      self.lims = [tmin,tmax]
      if self.boxh:
        self.boxh.remove()
      if self.direction is 'horizontal':
        self.boxh = pl.axvspan(tmin,tmax,facecolor='r',alpha=0.5)
	delta = tmaz-tmin
	axh.set_xlim([tmin-0.2*deltai, tmax+0.2*delta])
	if data is not None:
	  dat = data[1]

      if self.direction is 'vertical':
        self.boxh = pl.axhspan(tmin,tmax,facecolor='r',alpha=0.5)
	axh.set_ylim([tmin-0.2*delta, tmax+0.2*delta])
      fig.canvas.draw()
    
    def button_press_callback(self,event):
      if event.inaxes:
        if event.button == 3:
          self.finished = True
  roi = ROI(fig,axh,direction)
  selector = pl.matplotlib.widgets.SpanSelector(axh,roi.coo,direction)
  fig.canvas.mpl_connect('button_press_event', roi.button_press_callback)
  print "Select Span region of interest, finish with right click."
  while not roi.finished:
    pl.waitforbuttonpress()
  print "Span %s selected."%(roi.lims) 
  roi.boxh.remove()
  fig.canvas.draw()
  del selector
  return roi.lims

def getBins(lims,direction='horizontal',axh=None,fig=None):
  # TODO
  """Tool for selecting a span, functionality similar to ginput. Finish with right mouse button."""
  if not axh:
    axh = pl.gca()
  if not fig: fig = pl.gcf()
  class BINNING:
    def __init__(self,fig,axh,direction):
      self.fig = fig
      self.axh = axh
      self.lims = []
      self.boxh = []
      self.finished = False
      self.direction = direction

    def coo(self,tmin,tmax):
      self.lims = [tmin,tmax]
      if self.boxh:
        self.boxh.remove()
      if self.direction is 'horizontal':
        self.boxh = pl.axvspan(tmin,tmax,facecolor='r',alpha=0.5)
      if self.direction is 'vertical':
        self.boxh = pl.axhspan(tmin,tmax,facecolor='r',alpha=0.5)
      fig.canvas.draw()
    
    def button_press_callback(self,event):
      if event.inaxes:
        if event.button == 3:
          self.finished = True
  roi = ROI(fig,axh,direction)
  selector = pl.matplotlib.widgets.SpanSelector(axh,roi.coo,direction)
  fig.canvas.mpl_connect('button_press_event', roi.button_press_callback)
  print "Select Span region of interest, finish with right click."
  while not roi.finished:
    pl.waitforbuttonpress()
  
  roi.boxh.remove()
  fig.canvas.draw()
  del selector
  return roi.lims
def getCoordinate(direction='both',axh=None,fig=None):
  """Tool for selecting a coordinate, functionality similar to ginput for a single point. Finish with right mouse button."""
  if not axh:
    axh = pl.gca()
  if not fig: fig = pl.gcf()
  hor=False;ver=False
  if direction is 'horizontal' or 'hor' or 'both':
    hor=True
  if direction is 'vertical' or 'ver' or 'both':
    ver=True

  finished=False
  def button_press_callback(event):
    if event.inaxes:
      if event.button == 3:
        finished = True
  fig.canvas.mpl_connect('button_press_event', button_press_callback)
  print "Select a coordinate, finish with right click."
  linh = []
  while not finished:
    for tlinh in linh:
      tlinh.remove()
      linh = []
    pl.draw()
    pos = pl.ginput(1)[0]
    if hor:
      linh.append(pl.axvline(pos[0]))
    if ver:
      linh.append(pl.axhline(pos[1]))
    pl.draw()
    pl.waitforbuttonpress()

  
  fig.canvas.draw()
  return pos


def getRectangleCoordinates(axh=None,fig=None):
  """Tool for selecting a rectangle, functionality similar to ginput. Finish with right mouse button."""
  if not axh:
    axh = pl.gca()
  if not fig: fig = pl.gcf()
  class ROI:
    def __init__(self,fig,axh):
      self.fig = fig
      self.axh = axh
      self.lims = []
      self.boxh = []
      self.finished = False

    def coo(self,eclick,erelease):

      self.lims = [min([eclick.xdata,erelease.xdata]),
                   max([eclick.xdata,erelease.xdata]),
                   min([eclick.ydata,erelease.ydata]),
                   max([eclick.ydata,erelease.ydata])]
                   
      if self.boxh:
        self.boxh.remove()
      ptch = pl.Rectangle([self.lims[0],self.lims[2]],self.lims[1]-self.lims[0],self.lims[3]-self.lims[2],facecolor='r',alpha=0.5,ec='k')
      self.boxh = self.axh.add_patch(ptch)
      fig.canvas.draw()
    
    def button_press_callback(self,event):
      if event.inaxes:
        if event.button == 3:
          self.finished = True
  roi = ROI(fig,axh)
  selector = pl.matplotlib.widgets.RectangleSelector(axh,roi.coo)
  fig.canvas.mpl_connect('button_press_event', roi.button_press_callback)
  print "Select rectangular region of interest, finish with right click."
  while not roi.finished:
    pl.waitforbuttonpress()
  
  roi.boxh.remove()
  del selector
  axh.patches[-1].remove()
  fig.canvas.draw()
  return roi.lims


def hist2dSmart(x,y,fac=400.,**kwargs):
  h,xe,ye = histogram2dSmart(x,y,fac,**kwargs)
  imagesc(xe,ye,h)


def histSmart(x,linespec='k',label='histSmart',fac=20.,**kwargs):
  y,x = histogramSmart(x,fac,**kwargs)
  pl.step(np.hstack([x,x[-1]]),np.hstack([0,y,0]),linespec,label=label)

def imagesc(*args,**kwargs):
  slider = kwargs.get('slider',False)
  handle = kwargs.get('handle',[])
  axes   = kwargs.get('axes',pl)
  if ("axes" in kwargs): del kwargs["axes"]
  if slider and handle==[]:
    pl.clf()
  interpolation = 'nearest'
  if len(args)==1:
    I = args[0]
    (ny,nx)=I.shape
    x = np.arange(nx)
    y = np.arange(ny)
  elif len(args)==3:
    x = args[0]
    y = args[1]
    I = args[2]
    ny,nx = np.shape(I)
    if len(x)==nx+1:
      print "imagesc: x-vector is one element longer than corresponding intensity matrix. Is interpreted as bin edges."
      x = x[:-1]+np.mean(np.diff(x))
    elif len(x)==nx-1:
      print "imagesc: x-vector is one element shorter than corresponding intensity matrix. Is interpreted as bin edges with overflow on both sides."
      x = x-np.mean(np.diff(x))
    elif not len(x)==nx:
      print "Warning! x vector does not fit the side length of your matrix, possibly transposed!"
    if len(y)==ny+1:
      print "imagesc: y-vector is one element longer than corresponding intensity matrix. Is interpreted as bin edges."
      y = y[:-1]+np.mean(np.diff(y))
    elif len(y)==ny-1:
      print "imagesc: y-vector is one element shorter than corresponding intensity matrix. Is interpreted as bin edges with overflow on both sides."
      y = y-np.mean(np.diff(y))
    elif not len(y)==ny:
      print "Warning! x vector does not fit the side length of your matrix, possibly transposed!"
      print "Warning! x vector does not fit the side length of your matrix, possibly transposed!"
  dx = float(max(x)-min(x))/nx
  dy = float(max(y)-min(y))/ny
  xmn = min(x)-dx/2
  xmx = max(x)+dx/2
  ymn = min(y)-dy/2
  ymx = max(y)+dy/2
  def _format_coord(x, y):
    icol = int( (x-xmn)/(dx)+0.5 ) 
    irow = int( (y-ymn)/(dy)+0.5 )
    try:
      z = I[irow,icol]
      return 'x=%1.4f, y=%1.4f, z=%1.4f'%(x, y, z)
    except IndexError:
      return 'x=%1.4f, y=%1.4f'%(x, y)
  h = axes.imshow(I,interpolation=interpolation,
              extent=[xmn,xmx,ymn,ymx],origin='bottom',**kwargs)
  h.axes.format_coord=_format_coord
  axes.axis('tight')
  
  slider=False
  if slider:
    from matplotlib.widgets import Slider, Button, RadioButtons
    fig = pl.gcf()
    fig.subplots_adjust(left=0.25, bottom=0.25)
    fig.colorbar(h)
    axcolor = 'lightgoldenrodyellow'
    prc = np.percentile(I,range(0,101,10))
    axmin = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    axmax = fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    for n in range(1,11):
      axmin.axvline(prc[n],color='k')
      axmax.axvline(prc[n],color='k')
    smin = Slider(axmin, 'Min', prc[0], prc[-1], valinit=prc[3])
    smax = Slider(axmax, 'Max', prc[0], prc[-1], valinit=prc[-4])
    def update(val):
      h.set_clim([smin.val,smax.val])
      fig.canvas.draw()
    smin.on_changed(update)
    smax.on_changed(update)
  return h

def imagesc2(*args,**kwargs):
  slider  = kwargs.get('slider',False)
  handle = kwargs.get('handle',[])
  if slider and handle==[]:
    pl.clf()
  interpolation = 'nearest'
  if len(args)==1:
    h = pl.imshow(args[0],interpolation=interpolation)
    I = args[0]
  elif len(args)==3:
    x = args[0]
    y = args[1]
    I = args[2]
    ny,nx = np.shape(I)
    if len(x)==nx+1:
      print "imagesc: x-vector is one element longer than corresponding intensity matrix. Is interpreted as bin edges."
      x = x[:-1]+np.mean(np.diff(x))
    elif len(x)==nx-1:
      print "imagesc: x-vector is one element shorter than corresponding intensity matrix. Is interpreted as bin edges with overflow on both sides."
      x = x-np.mean(np.diff(x))
    elif not len(x)==nx:
      print "Warning! x vector does not fit the side length of your matrix, possibly transposed!"
    if len(y)==ny+1:
      print "imagesc: y-vector is one element longer than corresponding intensity matrix. Is interpreted as bin edges."
      y = y[:-1]+np.mean(np.diff(y))
    elif len(y)==ny-1:
      print "imagesc: y-vector is one element shorter than corresponding intensity matrix. Is interpreted as bin edges with overflow on both sides."
      y = y-np.mean(np.diff(y))
    elif not len(y)==ny:
      print "Warning! x vector does not fit the side length of your matrix, possibly transposed!"
      print "Warning! x vector does not fit the side length of your matrix, possibly transposed!"
    dx = float(max(x)-min(x))
    dy = float(max(y)-min(y))
    xmn = min(x)-dx/(nx-1)/2
    xmx = max(x)+dx/(nx-1)/2
    ymn = min(y)-dy/(ny-1)/2
    ymx = max(y)+dy/(ny-1)/2
    h = pl.imshow(I,interpolation=interpolation,origin='bottom',**kwargs)
    (loc,lab) = pl.xticks()
    idx = (loc>=0) & (loc<nx)
    loc = loc[idx]
    pl.xticks(loc,[x[p] for p in loc])
    pl.axis('tight')

  if slider:
    from matplotlib.widgets import Slider, Button, RadioButtons
    fig = pl.gcf()
    fig.subplots_adjust(left=0.25, bottom=0.25)
    fig.colorbar(h)
    axcolor = 'lightgoldenrodyellow'
    prc = np.percentile(I,range(0,101,10))
    axmin = fig.add_axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
    axmax   = fig.add_axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
    for n in range(1,11):
      axmin.axvline(prc[n],color='k')
      axmax.axvline(prc[n],color='k')
    smin = Slider(axmin, 'Min', prc[0], prc[-1], valinit=prc[3])
    smax = Slider(axmax, 'Max', prc[0], prc[-1], valinit=prc[-4])
    def updatemax(val):
      if smax.val<smin.val:
	smax.set_val(smin.val)
      h.set_clim([smin.val,smax.val])
      fig.canvas.draw()
    def updatemin(val):
      if smax.val<smin.val:
	smin.set_val(smax.val)
      h.set_clim([smin.val,smax.val])
      fig.canvas.draw()
    smin.on_changed(updatemin)
    smax.on_changed(updatemax)

  return h

def adjustYlabels(axes):
  pass


def subplots(*args,**kwargs):
  figname = kwargs.pop('figname',None)
  hspacemin = kwargs.pop('hspacemin',0.05)
  hspace = kwargs.pop('hspace',None)
  top = kwargs.pop('top',None)
  bottom = kwargs.pop('bottom',None)
  left = kwargs.pop('left',None)
  right = kwargs.pop('right',None)
  wspace = kwargs.pop('wspace',None)
  sharex = kwargs.get('sharex',False)
  if figname is not None:
    kwargs['num']=nfigure(figname).number
  fig,ax = plt.subplots(*args,**kwargs)
  plt.subplots_adjust(hspace=hspace,bottom=bottom,left=left,right=right,top=top,wspace=wspace)
  if (fig.subplotpars.hspace < hspacemin) and sharex:
    print 'YES'
    #for tax in ax[:-1]:
      #tax.set_xticklabels([])
      #xticklabels += tax.get_xticklabels()
    #plt.setp(xticklabels, visible=False)
    nbins = len(ax[0].get_xticklabels()) # added 
    for tax in ax[1:]:
      tax.yaxis.set_major_locator(MaxNLocator(nbins=nbins, prune='upper')) # added 
  return fig,ax


def addPngMetadata(fname,dictionary):
	""" This function is meant to store information about the analysis on a 
	png image created with matplotlib;
	It takes a filename (of a previously saved png image) and a dictionary;
	inspired by Nick Galbreath google code: pngaddcomment
	"""
	outtemp = "addPngMetadata.png"
	try:
		from PIL import Image
		from PIL import PngImagePlugin
	except ImportError:
		print "PIL module not found, can't add PngMetadata"
	img = Image.open(fname)
	meta = PngImagePlugin.PngInfo()
	for k,v in dictionary.iteritems():
		meta.add_text(k,str(v))
#	for k,v in img.info.iteritems():
#		meta.add_text(k,str(v))
	img.save(outtemp,"PNG", pnginfo=meta)
# if previous steps works, writing directly to fname could corrupt file
	os.rename(outtemp,fname);
	 
def getPngMetadata(fname):
	""" This function is meant to store information about the analysis on a 
	png image created with matplotlib;
	It takes a filename (of a previously saved png image) and a dictionary;
	inspired by Nick Galbreath google code: pngaddcomment
	"""
	try:
		from PIL import Image
		from PIL import PngImagePlugin
	except ImportError:
		print "PIL module not found, can't add PngMetadata"
		return None
	img = Image.open(fname)
	return img.info


from matplotlib.colors import LinearSegmentedColormap

class nlcmap(LinearSegmentedColormap):
    """A nonlinear colormap"""
    name = 'nlcmap'
    def __init__(self, cmap, levels):
        self.cmap = cmap
        # @MRR: Need to add N for backend
        self.N = cmap.N
        self.monochrome = self.cmap.monochrome
        self.levels = np.asarray(levels, dtype='float64')
        self._x = self.levels / self.levels.max()
        self._y = np.linspace(0.0, 1.0, len(self.levels))
    #@MRR Need to add **kw for 'bytes'
    def __call__(self, xi, alpha=1.0, **kw):
        """docstring for fname"""
        # @MRR: Appears broken?
        # It appears something's wrong with the
        # dimensionality of a calculation intermediate
        #yi = stineman_interp(xi, self._x, self._y)
        yi = np.interp(xi, self._x, self._y)
        return self.cmap(yi, alpha)
 
 
def cmap_smooth(I=None,axh=None,Nlevels=256,cmap_lin=None):
  if cmap_lin==None:
    cmap_lin = pl.cm.jet
  if I == None:
    if not axh:
      axh = pl.gca()
    ihandles = axh.findobj(matplotlib.image.AxesImage)
    ih = ihandles[-1]
    I = ih.get_array()
  levels = np.percentile(I.ravel(),list(np.linspace(0,100,Nlevels)))
  cmap_nonlin = nlcmap(cmap_lin,levels)
  pl.set_cmap(cmap_nonlin)


greyscale = [
            " ",
            " ",
            ".,-",
            "_ivc=!/|\\~",
            "gjez2]/(YL)t[+T7Vf",
            "mdK4ZGbNDXY5P*Q",
            "W8KMA",
            "#%$"
            ]
from bisect import bisect
from random import randint

def hist_asciicontrast(x,bins=50,range=None,disprange=True):
  h,edges = np.histogram(x,bins=bins,range=range)
  if np.sum(h)==0:
    bounds = np.linspace(min(h),1,len(greyscale)-1)
  else:
    bounds = np.linspace(min(h),max(h),len(greyscale)-1)
  hstr = ''
  
  for bin in h:
    syms = greyscale[bisect(bounds,bin)]
    hstr+=syms[randint(0,len(syms)-1)]

  if disprange:
    hstr = '{:>10}'.format('%0.5g' %(edges[0])) + hstr + '{:>10}'.format('%0.5g' %(edges[-1]))

  return hstr  
  
class hist_ascii(object):
    """
    Ascii histogram
    """
    def __init__(self, data, bins=50,percRange=None,range=None):
        """
        Class constructor
        
        :Parameters:
            - `data`: array like object
        """
	if not percRange==None:
	  range = np.percentile(data,percRange)

        self.data = data
        self.bins = bins
        self.h = np.histogram(self.data, bins=self.bins, range=range)
    def horizontal(self, height=4, character ='|'):
        """Returns a multiline string containing a
        a horizontal histogram representation of self.data
        :Parameters:
            - `height`: Height of the histogram in characters
            - `character`: Character to use
        >>> d = normal(size=1000)
        >>> h = Histogram(d,bins=25)
        >>> print h.horizontal(5,'|')
        106            |||
                      |||||
                      |||||||
                    ||||||||||
                   |||||||||||||
        -3.42                         3.09
        """
        his = """"""
        bars = 1.*self.h[0]/np.max(self.h[0])*height
        formnum = lambda(num): '{:<9}'.format('%0.4g' %(num))
        for l in reversed(range(1,height+1)):      
            line = ""
            if l == height:
                line = formnum(np.max(self.h[0])) + ' ' #histogram top count
            else:
                line = ' '*(9+1) #add leading spaces
            for c in bars:
                if c >= np.ceil(l):
                    line += character
                else:
                    line += ' '
            line +='\n'
            his += line
        his += formnum(self.h[1][0]) + ' '*(self.bins) + formnum(self.h[1][-1]) + '\n'
        return his
    def vertical(self,height=20, character ='|'):
        """
        Returns a Multi-line string containing a
        a vertical histogram representation of self.data
        :Parameters:
            - `height`: Height of the histogram in characters
            - `character`: Character to use
        >>> d = normal(size=1000)
        >>> Histogram(d,bins=10)
        >>> print h.vertical(15,'*')
                              236
        -3.42:
        -2.78:
        -2.14: ***
        -1.51: *********
        -0.87: *************
        -0.23: ***************
        0.41 : ***********
        1.04 : ********
        1.68 : *
        2.32 :
        """
        his = """"""
        xl = ['%.2f'%n for n in self.h[1]]
        lxl = [len(l) for l in xl]
        bars = self.h[0]/max(self.h[0])*height
        his += ' '*(np.max(bars)+2+np.max(lxl))+'%s\n'%np.max(self.h[0])
        for i,c in enumerate(bars):
            line = xl[i] +' '*(np.max(lxl)-lxl[i])+': '+ character*c+'\n'
            his += line
        return his
            
 

class Roipoly(Widget):
  def __init__(self, ax, xy=None, callback=None, useblit=True):
    self.axes = ax
    self.figure = ax.figure
    self.canvas = self.figure.canvas
    self.useblit = useblit
    if useblit:
	self.background = self.canvas.copy_from_bbox(self.axes.bbox)
    if not xy==None:
      x, y = xy
      self.verts = [(x,y)]
    else:
      self.verts = []

    self.callback = callback
    self.line=None
    self.dragging = False
    self.open_path = False
    self.done = False
    self.cids = []

    self.cids.append(self.canvas.mpl_connect('button_press_event', self.onclick))
    self.cids.append(self.canvas.mpl_connect('button_release_event', self.onrelease))
    self.cids.append(self.canvas.mpl_connect('motion_notify_event', self.onmove))

  def onrelease(self, event):
    if self.dragging:
      self.verts.append((event.xdata, event.ydata))
      self.dragging = False
      #if len(self.verts)>2:
	#self.callback(self.verts)
#

  def onclick(self, event):
    if event.inaxes != self.axes: return
    self.clicktime = time.time()
    if event.button==3:
      if self.open_path:
	self.open_path = False
        self.verts.append((event.xdata, event.ydata))
	plotdat = copy(self.verts)
	plotdat.append(plotdat[0])
	self.line.set_data(zip(*plotdat))
        self.axes.lines.remove(self.line_close)
	self.line_close=None
      else:
	
	for cid in self.cids:
	  self.canvas.mpl_disconnect(cid)
        self.axes.lines.remove(self.line)
	self.line = None
	self.done = True

	


      #self.verts = None
      #self.axes.lines.remove(self.line)
      #self.verts = None
      #for cid in self.cids:
	  #self.canvas.mpl_disconnect(cid)
      #return 
    if event.button==1:
      if not self.open_path:
        self.verts = []

	self.line=None
	#self.dragging = False
	#self.open_path = False

      self.verts.append((event.xdata, event.ydata))
      if self.line==None:
	self.verts.append((event.xdata, event.ydata))
	self.open_path = True
	x,y = zip(*self.verts)
	self.line = Line2D([x], [y], linestyle='-', color='black', lw=1.3)
	self.axes.add_line(self.line)
	self.line_close = Line2D([x], [y], linestyle=':', color='black', lw=1.3)
	self.axes.add_line(self.line_close)

    if self.useblit:
	self.canvas.restore_region(self.background)
	if not self.line==None:
	  self.axes.draw_artist(self.line)
	if not self.line_close==None:
	  self.axes.draw_artist(self.line_close)
	self.canvas.blit(self.axes.bbox)
    else:
	self.canvas.draw_idle()

  def onmove(self, event):
    if not self.open_path: return
    if self.verts==[]: return
    if event.inaxes != self.axes: return
    if event.button==1 and (time.time()-self.clicktime)>.1:
      self.dragging = True
      self.verts.append((event.xdata, event.ydata))
      plotdat = copy(self.verts)
      closedat = [plotdat[-1],plotdat[0]]
      self.line.set_data(zip(*plotdat))
      self.line_close.set_data(zip(*closedat))
    elif event.button==None:
      plotdat = copy(self.verts)
      plotdat.append((event.xdata, event.ydata))
      closedat = [plotdat[-1],plotdat[0]]
      self.line.set_data(zip(*plotdat))
      self.line_close.set_data(zip(*closedat))

    if self.useblit:
	self.canvas.restore_region(self.background)
	self.axes.draw_artist(self.line)
	self.axes.draw_artist(self.line_close)
	self.canvas.blit(self.axes.bbox)
    else:
	self.canvas.draw_idle()

def roipoly(ax=None):
  """ use polygonmask to get the pixel mask"""
  if ax==None:
    ax = plt.gca()
  r = Roipoly(ax)
  while not r.done:
    plt.waitforbuttonpress()
  return r.verts

def polygonmask(pgon,x,y):
  xshp = np.shape(x)
  yshp = np.shape(y)
  if not xshp==yshp:
    print "Shapes dont match"
    return
  if LooseVersion(matplotlib.__version__)<'1.2':
    mask = points_inside_poly(zip(x.ravel(),y.ravel()),pgon)
  else:
    tpath = mPath(pgon) 
    mask = tpath.contains_points(zip(x.ravel(),y.ravel()))
  

  mask = np.reshape(mask,xshp)
  return mask

def colormaps():
  from matplotlib.colors import LinearSegmentedColormap

  """

  Example: suppose you want red to increase from 0 to 1 over the bottom
  half, green to do the same over the middle half, and blue over the top
  half.  Then you would use:

  cdict = {'red':   ((0.0,  0.0, 0.0),
		     (0.5,  1.0, 1.0),
		     (1.0,  1.0, 1.0)),

	   'green': ((0.0,  0.0, 0.0),
		     (0.25, 0.0, 0.0),
		     (0.75, 1.0, 1.0),
		     (1.0,  1.0, 1.0)),

	   'blue':  ((0.0,  0.0, 0.0),
		     (0.5,  0.0, 0.0),
		     (1.0,  1.0, 1.0))}

  If, as in this example, there are no discontinuities in the r, g, and b
  components, then it is quite simple: the second and third element of
  each tuple, above, is the same--call it "y".  The first element ("x")
  defines interpolation intervals over the full range of 0 to 1, and it
  must span that whole range.  In other words, the values of x divide the
  0-to-1 range into a set of segments, and y gives the end-point color
  values for each segment.

  Now consider the green. cdict['green'] is saying that for
  0 <= x <= 0.25, y is zero; no green.
  0.25 < x <= 0.75, y varies linearly from 0 to 1.
  x > 0.75, y remains at 1, full green.

  If there are discontinuities, then it is a little more complicated.
  Label the 3 elements in each row in the cdict entry for a given color as
  (x, y0, y1).  Then for values of x between x[i] and x[i+1] the color
  value is interpolated between y1[i] and y0[i+1].

  Going back to the cookbook example, look at cdict['red']; because y0 !=
  y1, it is saying that for x from 0 to 0.5, red increases from 0 to 1,
  but then it jumps down, so that for x from 0.5 to 1, red increases from
  0.7 to 1.  Green ramps from 0 to 1 as x goes from 0 to 0.5, then jumps
  back to 0, and ramps back to 1 as x goes from 0.5 to 1.

  row i:   x  y0  y1
		  /
		 /
  row i+1: x  y0  y1

  Above is an attempt to show that for x in the range x[i] to x[i+1], the
  interpolation is between y1[i] and y0[i+1].  So, y0[0] and y1[-1] are
  never used.

  """



  cdict1 = {'red':   ((0.0, 0.0, 0.0),
		     (0.5, 0.0, 0.1),
		     (1.0, 1.0, 1.0)),

	   'green': ((0.0, 0.0, 0.0),
		     (1.0, 0.0, 0.0)),

	   'blue':  ((0.0, 0.0, 1.0),
		     (0.5, 0.1, 0.0),
		     (1.0, 0.0, 0.0))
	  }

  cdict2 = {'red':   ((0.0, 0.0, 0.0),
		     (0.5, 0.0, 1.0),
		     (1.0, 0.1, 1.0)),

	   'green': ((0.0, 0.0, 0.0),
		     (1.0, 0.0, 0.0)),

	   'blue':  ((0.0, 0.0, 0.1),
		     (0.5, 1.0, 0.0),
		     (1.0, 0.0, 0.0))
	  }

  cdict3 = {'red':  ((0.0, 0.0, 0.0),
		     (0.25,0.0, 0.0),
		     (0.5, 0.8, 1.0),
		     (0.75,1.0, 1.0),
		     (1.0, 0.4, 1.0)),

	   'green': ((0.0, 0.0, 0.0),
		     (0.25,0.0, 0.0),
		     (0.5, 0.9, 0.9),
		     (0.75,0.0, 0.0),
		     (1.0, 0.0, 0.0)),

	   'blue':  ((0.0, 0.0, 0.4),
		     (0.25,1.0, 1.0),
		     (0.5, 1.0, 0.8),
		     (0.75,0.0, 0.0),
		     (1.0, 0.0, 0.0))
	  }

  # Make a modified version of cdict3 with some transparency
  # in the middle of the range.
  cdict4 = cdict3.copy()
  cdict4['alpha'] = ((0.0, 1.0, 1.0),
		  #   (0.25,1.0, 1.0),
		     (0.5, 0.3, 0.3),
		  #   (0.75,1.0, 1.0),
		     (1.0, 1.0, 1.0))


  # Now we will use this example to illustrate 3 ways of
  # handling custom colormaps.
  # First, the most direct and explicit:

  blue_red1 = LinearSegmentedColormap('BlueRed1', cdict1)

  # Second, create the map explicitly and register it.
  # Like the first method, this method works with any kind
  # of Colormap, not just
  # a LinearSegmentedColormap:

  blue_red2 = LinearSegmentedColormap('BlueRed2', cdict2)
  plt.register_cmap(cmap=blue_red2)

  # Third, for LinearSegmentedColormap only,
  # leave everything to register_cmap:

  plt.register_cmap(name='BlueRed3', data=cdict3) # optional lut kwarg
  plt.register_cmap(name='BlueRedAlpha', data=cdict4)


def showRealTimeDevelopment(*args,**kwargs):
  figname = kwargs.get('figname','realtime_development')
  resolution = kwargs.get('resolution',1e-9)
  plotIt = kwargs.get('plotIt',True)
  N = len(args)
  if plotIt:
    fig = nfigure(figname)
    fig,ax = plt.subplots(N,1,sharex=True,num=fig.number)
  opdata = []
  for tax,var in zip(ax,args):
    t = np.hstack(var.getTime(resolution=resolution,asTime=True))
    if plotIt:
      tax.plot(t,var.R,'.k',ms=1)
    opdata.append([t,var.R])
  return opdata



