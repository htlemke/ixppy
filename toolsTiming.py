from scipy import optimize,special
import numpy as np
import utilities
import sys
sqrt2=np.sqrt(2.)

def TTfuncFit(x,x0,a,sig,b0,b1):
  sig=abs(sig)
  step = a*(special.erf((x-x0)/sig/sqrt2)+1)/2
  bkg  = (b0+b1*x)
  return step+bkg

def findDerPeak(x,der,excludePoints=50,use="max"):
    if (use=="max"):
      arg = der[excludePoints:-excludePoints].argmax(); # exclude points at extreme
    else:
      arg = der[excludePoints:-excludePoints].argmin(); # exclude points at extreme
    return x[arg+excludePoints]


def findStepWithPoly(x,data,kind="stepUp",excludePoints=100,order=20,fitrange=100):
	""" Look for a step in the data
	    Data is 1D array
	    The 'kind' keywords should be either 'stepUp' or 'stepDown'
	    the 'excludePoints' keyword is used to limit the search 'excludePoints' 
	    away from the extremes
	    The data are fit with a polynomial of order 'order', then the maximum
	    (or minimum) of derivative is located.
	    After this first attempt, the position is refined by a second polyfit
			in a range [-fitrange,+fitrange] around the first guess
	"""
	if (kind == "stepUp"):
		use = "max"
	else:
		use = "min"

	poly    = np.polyfit(x,data,order)
	polyder = np.polyder(poly)
	x_poly1 = findDerPeak(x,np.polyval(polyder,x),use=use,excludePoints=excludePoints)
	idx = ( x>(x_poly1-fitrange) ) & (x<(x_poly1+fitrange) )
	poly    = np.polyfit(x[idx],data[idx],order)
	polyder = np.polyder(poly)
	x_poly2 = findDerPeak(x[idx],np.polyval(polyder,x[idx]),use=use,excludePoints=10)
	return x_poly2

def findStepWithErfFit(x,data,kind="stepUp",excludePoints=100,order=20,fitrange=50):
	""" Look for a step in the data
	    Data can be a 1D array or a 2D ones, in the latter case the 'axis' index
	    if used as different shot index
	    The 'kind' keywords should be either 'stepUp' or 'stepDown'
	    the 'excludePoints' keyword is used to limit the search 'excludePoints' 
	    away from the extremes
	"""
	x_poly = findStepWithPoly(x,data,kind=kind,excludePoints=excludePoints,order=order,fitrange=fitrange)
	idx = ( x>(x_poly-fitrange) ) & (x<(x_poly+fitrange) )
	xfit = x[idx]; y = data[idx]
	# estimate errors by high order polinomial fit
	p = np.polyfit(xfit,y,30)
	err = y-np.polyval(p,xfit)
	err = np.std(err)
	# autoguess parameters
	sig = fitrange/6.
	meanLeft = np.mean(y[:10])
	meanRight= np.mean(y[-10:])
	a   = meanRight-meanLeft
	b0  = meanLeft
	b1  = 0.001
	fitp=optimize.curve_fit(TTfuncFit,xfit,y,p0=(x_poly,a,sig,b0,b1),\
	maxfev=10000,ftol=1e-3, sigma=err)
	(x0,a,sig,b0,b1) = fitp[0]
	try:
		(ex0,ea,esig,eb0,eb1) = np.sqrt( np.diag( fitp[1] ) )
	except:
		(ex0,ea,esig,eb0,eb1) = (0,0,0,0,0)
	if (x0>xfit.max()) or (x0<xfit.min()): x0 = xfit.mean()
	yfit = TTfuncFit(xfit,x0,a,sig,b0,b1)
	return x0,a,sig,xfit,yfit,ex0,ea,esig

def findStepImage(img,img_bkg=None,roi=None,roiref=None,roinorm=None,run=None,axis=0,use="max"):
	""" find the step in the roi of the image IM obtained as:
				IM = img/roinorm.mean(img) -
						 img_bkg/roinorm.mean(img_bkg), if img_bkg is not None
			or
				IM = roi.select(img)/roinorm.mean( roi.select(img) ) -
					 - roiref.select(img)/roinorm.mean( roiref.select(img) )
						 if img_bkg is None or roiref is not None
			In the latter case the size of roi and roiref has to be the same
			If roiref is given, roinorm is defined in the roi frame
			ROIs have to be passed or alternatively the run  """
	if (run is not None): 
		(roi,roiref,roinorm) = runToROIs(run)
	elif (roi is None):
		print "Either run or rois has to be defined, exiting"
		sys.exit(1)
	# Calculate image and background, normalizing if needed
	# note that the normalization is done differently if roiref is given or not
	if   ( (img_bkg is None) or (roiref is not None) ):
		(img,img_bkg) = (roi.select(img),roiref.select(img))
		if (roinorm is not None):
			img     /= roinorm(img)
			img_bkg /= roinorm(img_bkg)
	else:
		if (roinorm is not None):
			n     = roinorm(img)
			n_bkg = roinorm(img_bkg)
		else:
			n=n_bkg=1
 		(img,img_bkg) = (roi.select(img),roi.select(img_bkg))
		img /= n; img_bkg/=n_bkg;
	diff_img   = img-img_bkg
	diff_curve = img_diff.mean(axis=axis)
	x = np.arange(roi.cmin,roi.cmax)
	return findStepCurve(c,diff_curve,use=use)
