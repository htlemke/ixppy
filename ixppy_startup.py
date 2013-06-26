import os,sys
#sys.stdout.flush()
print '\n\n'
print '.__                             '
print '|__|__  _________ ______ ___.__.'
print '|  \  \/  /\____ \\\\____ <   |  |'
print '|  |>    < |  |_> >  |_> >___  |'
print '|__/__/\_ \|   __/|   __// ____|'
print '         \/|__|   |__|   \/     '
print 'Interactive X-ray Pulse Python Environemt\n\n'

import ixppy
import h5py
parentdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,parentdir) 
from tools import *
get_ipython().run_line_magic('load_ext','autoreload')
get_ipython().run_line_magic('autoreload','2')
h5py.enable_ipython_completer()
