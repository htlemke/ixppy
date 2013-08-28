import numpy as np
from toolsConstsAndConv import *
from toolsVarious import *

def getLODCMdelay(E,oldE,ds=.6,ID='Si',hkl=(1,1,1)):
    theold = BraggAngle(ID,hkl,E=oldE)
    thenew = BraggAngle(ID,hkl,E=E)
    pathold = ds/sind(2*theold)-ds/tand(2*theold)
    pathnew = ds/sind(2*thenew)-ds/tand(2*thenew)
    delay = (pathnew-pathold)/c_light()
    return delay

