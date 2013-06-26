from pylab import *
import ixppy
import sys
import os,sys
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,parentdir) 
from tools import *


def read_data():
  d = ixppy.dataset('/reg/g/xpp/data/example_data/diode/time_tool/hdf5/xppi0412-r0113.h5',['timeTool','ipm2','ipm3','diodeU','opal0','eventCode'])
  
  # scanning motor
  print d.scanMot

  # number of scansteps
  print d.noofsteps
  
  return d

def monitor(d=None):
  if d==None: d = read_data()
  # look at ipm data
  # stack all data together
  nfigure('IPMs')
  ipm2 = hstack(d.ipm2.sum)
  hist(ipm2,500)

def raw_data(d=None):
  if d==None: d = read_data()
  # plot average of data
  nfigure('average data')
  I = [ mean(i) for i in d.diodeU.channel1]
  plot(d.scanVec,I)

def filter(d=None):
  if d==None: d = read_data()
  # filter function, graphical input if no limits passed
  ixppy.parameterFilt(d.ipm2.sum,d)
  return d

def plot_normalized(d=None):
  if d==None: d = read_data()
  #plot normalized data
  nfigure('average data')
  I = [ mean(i)/mean(i0) for i,i0 in zip(d.diodeU.channel1,d.ipm2.sum)]
  plot(d.scanVec,I)

def show_opal(d=None):
  if d==None: d = read_data()
  I = d.opal0.rdStepData(0,range(20))
  print shape(I)
  Imean = mean(I,axis=0)
  nfigure('Opal image')
  imagesc(Imean)

  # chosen vertical lines
  roilims = [540,590]
  # Reduce dimension for working with the dataset
  I = squeeze(Imean)
  #plot selection
  nfigure('opal ROI')
  subplot(211)
  imshow(I[roilims[0]:roilims[1],:])
  axis('normal')
  axis('tight')
  # plot profile
  subplot(212)
  plot(mean(I[roilims[0]:roilims[1],:],0))

def timeTool_calibration():
  d = ixppy.dataset(('xppi0412',115),['timeTool'])
  nfigure('step Timetool hist')
  for t in d.timeTool.fltpos:
    clf()
    hist(t)

def chunks(d=None):
  if d==None: d = read_data()
  # generating chunks for larger scale data processing
  # Using Nmax argument for# generating chunks for larger scale data processing
  # Using Nmax argument for testing of code, run for 10 images per chunk first to see if it works.
  chunks = d.opal0.chunks(Nmax=10)
  roilims = [540,590]
  roilims = [590,630]
  allprofiles = []
  for calibc in chunks:
    ccprofiles = []
    for memchunk in calibc:
      #print shape(memchunk.data)
      this_anadata = squeeze(mean(memchunk.data[:,roilims[0]:roilims[1],:],1))
      #print shape(this_anadata)
      ccprofiles.append(this_anadata)
    # stack chunks together and append them to the final dataset.
    allprofiles.append(vstack(ccprofiles))
  return allprofiles

def chunks_hdf5(d=None):
  if d==None: d = read_data()
  # generating chunks for larger scale data processing
  # Using Nmax argument for# generating chunks for larger scale data processing
  # Using Nmax argument for testing of code, run for 10 images per chunk first to see if it works.
  chunks = d.opal0.chunks(Nmax=10)
  roilims = [540,590]
  roilims = [590,630]
  allprofiles = []
  for calibc in chunks:
    ccprofiles = []
    for memchunk in calibc:
      #print shape(memchunk.data)
      this_anadata = squeeze(mean(memchunk.data[:,roilims[0]:roilims[1],:],1))
      #print shape(this_anadata)
      ccprofiles.append(this_anadata)
    # stack chunks together and append them to the final dataset.
    d.opal0.append_h5_datalist('profile', vstack(ccprofiles))
  return d

def timeTool_binning(d=None):
  if d==None: d = read_data()
  ttcalib = [5.3106e-18,-1.2164e-14,5.687e-12]

  # filtering for the dropped shots
  filt = []
  for step in d.eventCode.code162:
    filt.append(step==True)
  d.filter = filt
  ixppy.parameterFilt(d.ipm2.sum,d,lims=[.1,10])
  # getting all important data
  i0 = hstack(d.ipm2.sum)
  i = hstack(d.diodeU.channel1)
  tt = hstack([delay+polyval(ttcalib,fltpos) for delay,fltpos in zip(d.scanVec,d.timeTool.fltpos)])
  hist2DSmart(tt,i/i0)
  return i,i0,tt

