print '\n\n'
print '.__                             '
print '|__|__  _________ ______ ___.__.'
print '|  \  \/  /\____ \\\\____ <   |  |'
print '|  |>    < |  |_> >  |_> >___  |'
print '|__/__/\_ \|   __/|   __// ____|'
print '         \/|__|   |__|   \/     '
print 'Interactive X-ray Pulse Python Environment\n\n'
import numpy as np
import matplotlib.pyplot as plt
import ixppy
from tools import *
import os,sys
import h5py
try:
  h5py.enable_ipython_completer()
except:
  print "h5py tab completion will not work in your python environment."
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,parentdir) 
get_ipython().run_line_magic('load_ext','autoreload')
get_ipython().run_line_magic('autoreload','1')
