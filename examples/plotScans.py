import sys
import ixppy
import numpy as np
import pylab as p
from ixppy import kulletools as kt
import utils
from matplotlib.backends.backend_pdf import PdfPages
g_outdir = "../analysis_out/"
PdfPlots = PdfPages(g_outdir+'plotScans.pdf')
#p.interactive(1)



DOCAL = True
nonLinTra = utils.correctNonLin()
nonLinFluo= utils.correctNonLin()

def calcScan(run,correctnonlin=True):
	dark_limit = 0.01
	d=ixppy.dataset(('xpp52512',run),['ipm2','ipm3',"diodeU","diode3"])
	d.filtTimestamps()
	Ncal	 = len(d.ipm2.sum)
	mon2 = d.ipm2.__sum()
	mon3 = d.ipm3.__sum()
	xpos = d.ipm3.__xpos()
	fluo = d.diodeU.__channel0()
	tra	= d.diode3.__channel0()
	mon = mon2
	x = d.scanVec
	if (d.scanMot == "ccmE_vernier"):
		correctnonlin=False
	kt.nfigure("Data run %s" % run,figsize=(8,4))
	p.clf()
	p.subplots_adjust(wspace=0.4, hspace=0.6)

	if (correctnonlin):
		M = np.hstack(mon[-4:-1])
		D = np.hstack(fluo[-4:-1])
		idx = (D<1.7) & (M>0.01) & (D>0.01) & (np.abs(D/M -np.median(D/M)) < 0.1)
		idx = (D<1.7) & (M>0.01); # & (D>0.01) & (np.abs(D/M -np.median(D/M)) < 0.1)
		nonLinFluo.use(M[idx],D[idx])
		D = np.hstack(tra[1:4])
		idx = (D<1.7) & (M>0.01) & (D>0.01) & (np.abs(D/M -np.median(D/M)) < 0.1)
		idx = (D<1.7) & (M>0.01);# & (D>0.01) & (np.abs(D/M -np.median(D/M)) < 0.1)
		nonLinTra.use(M[idx],D[idx])

	outF=[]
	outFinfo=[]
	outT=[]
	outTinfo=[]
	for i in range(Ncal):
		m		= mon[i]
		f		= fluo[i]
		t		= tra[i]
		A = utils.myAveragingAndFilter()
		A.addFilter( m>dark_limit,"mon > %s" % (dark_limit))
		A.addFilter( f<1.7,"fluo<1.7")
#	F.add( np.abs(TTs-np.median(TTs))<2*kt.mad(TTs),"TTs 1 sigma")
#		print F.info()
		if (correctnonlin):
			if (i==0):
				p.subplot("222",title="Non lin correction")
				p.plot(m,f/m,"+",label="fluo before");
				p.plot(m,t/m,"+",label="tra  before")
			(m,f) = nonLinFluo.correct(m,f)
			(m,t) = nonLinTra.correct(m,t)
			if (i==0):
				p.subplot("222",title="Non lin correction")
				p.plot(m,f/m,"o",label="fluo after")
				p.plot(m,t/m,"o",label="tra  after")
				p.legend(bbox_to_anchor = (0.8,-0.2))
				p.xlabel("mon")
				p.ylabel("[fluo|tra]/mon")
		(vF,infoF) = A.calcDiodeNormalized(m,f)
		outF.append(vF);
		outFinfo.append(infoF)
		(vT,infoT) = A.calcDiodeNormalized(m,t)
		outT.append(vT);
		outTinfo.append(infoT)
	f=open(g_outdir+"plotScan_run%s_fluo.txt" % run,"w")
	f.write("%s " % d.scanMot)
	H=A.header()
	for h in H:
		f.write("%s " % h)
	f.write("info\n")
	for i in range(Ncal):
		f.write("%s " % x[i])
		for j in range(len(outF[i])):
			f.write("%s " % (outF[i][j]))
		f.write("%s\n" % outFinfo[i])
	outF=np.array(outF)[:,0]
	outT=np.array(outT)[:,0]
	f.close()
	sub1 = p.subplot("221",title="Fluo run %s" % run)
	print x.shape,outF.shape
	p.plot(x,outF)
	p.subplot("223",title="Tra run %s" % run,sharex=sub1)
	p.xlabel(d.scanMot)
	p.plot(x,outT,label="run %s" % run)
	p.savefig(g_outdir+"plotScan_run%s.pdf" % run,trasparent=True,papertype="a4")
	PdfPlots.savefig(trasparent=True,papertype="a4")

#format = 	
#crystal2 = (34,185,"c2 hscan",
#						35,185,"c2 
runs = (213,214,215,216,218,219)
runs = (40,41,164,165)
runs = range(40,90)
for r in runs:
	print r
	d=ixppy.dataset(('xpp52512',r),['ipm2','ipm3',"diodeU","diode3"])
	if (hasattr(d,"scanMot") and (d.scanMot == "ccmE_vernier")):
		print "|||",r
		calcScan(r)
PdfPlots.close()
