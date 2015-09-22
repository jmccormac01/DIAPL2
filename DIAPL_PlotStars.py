# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_PlotStars.py - 	a program to filter out the variable stars
#										
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 31/10/11 - 	code writen
# 31/10/11 -	code tested
# 30/11/11 -	added outut for phased files to plot with gnuplot
# 08/12/11 - 	added error scaling to SExtractor mean
#				added errors to plots and outputs etc
#	
# 	Test1: Try Tamuzing!
#
#	Result: NO REAL CHANGE FROM TAMUZING!
#
# 	Test2: Address the systematics problem try several things. 
#
#			HJD values from header match the star files exactly!
#
# 		1 - only use observations from >40 deg in elevation
# 		2 - exclude any frame with large FWHM
#
# 		write a program to batch process each set, making a list of FWHM 
# 		and those images outside the desired altitude range. 
#
# 		then filter all lightcurves against this list according to their hjd value
# 		output a new list of files that have been filtered
#
# 		then finally test rejecting extreme outliers before binning
#	
#		Result: NO NOTICABLE CHANGE REMOVING THESE POINTS!

from pylab import *
import pyfits as pf
import commands, sys
import os, os.path
import numpy as np
from numpy import *
from numpy.fft import *

###############################################
################# FUNCTIONS ###################
###############################################

""" Fast algorithm for spectral analysis of unevenly sampled data

The Lomb-Scargle method performs spectral analysis on unevenly sampled
data and is known to be a powerful way to find, and test the 
significance of, weak periodic signals. The method has previously been
thought to be 'slow', requiring of order 10(2)N(2) operations to analyze
N data points. We show that Fast Fourier Transforms (FFTs) can be used
in a novel way to make the computation of order 10(2)N log N. Despite
its use of the FFT, the algorithm is in no way equivalent to 
conventional FFT periodogram analysis.

Keywords:
  DATA SAMPLING, FAST FOURIER TRANSFORMATIONS, 
  SPECTRUM ANALYSIS, SIGNAL  PROCESSING

Example:
  > import numpy
  > import lomb
  > x = numpy.arange(10)
  > y = numpy.sin(x)
  > fx,fy, nout, jmax, prob = lomb.fasper(x,y, 6., 6.)

Reference: 
  Press, W. H. & Rybicki, G. B. 1989
  ApJ vol. 338, p. 277-280.
  Fast algorithm for spectral analysis of unevenly sampled data
  bib code: 1989ApJ...338..277P

"""

def __spread__(y, yy, n, x, m):
  """
  Given an array yy(0:n-1), extirpolate (spread) a value y into
  m actual array elements that best approximate the "fictional"
  (i.e., possible noninteger) array element number x. The weights
  used are coefficients of the Lagrange interpolating polynomial
  Arguments:
    y : 
    yy : 
    n : 
    x : 
    m : 
  Returns:
    
  """
  nfac=[0,1,1,2,6,24,120,720,5040,40320,362880]
  if m > 10. :
    print 'factorial table too small in spread'
    return

  ix=long(x)
  if x == float(ix): 
    yy[ix]=yy[ix]+y
  else:
    ilo = long(x-0.5*float(m)+1.0)
    ilo = min( max( ilo , 1 ), n-m+1 ) 
    ihi = ilo+m-1
    nden = nfac[m]
    fac=x-ilo
    for j in range(ilo+1,ihi+1): fac = fac*(x-j)
    yy[ihi] = yy[ihi] + y*fac/(nden*(x-ihi))
    for j in range(ihi-1,ilo-1,-1):
      nden=(nden/(j+1-ilo))*(j-ihi)
      yy[j] = yy[j] + y*fac/(nden*(x-j))

def fasper(x,y,ofac,hifac, MACC=4):
  """ function fasper
    Given abscissas x (which need not be equally spaced) and ordinates
    y, and given a desired oversampling factor ofac (a typical value
    being 4 or larger). this routine creates an array wk1 with a
    sequence of nout increasing frequencies (not angular frequencies)
    up to hifac times the "average" Nyquist frequency, and creates
    an array wk2 with the values of the Lomb normalized periodogram at
    those frequencies. The arrays x and y are not altered. This
    routine also returns jmax such that wk2(jmax) is the maximum
    element in wk2, and prob, an estimate of the significance of that
    maximum against the hypothesis of random noise. A small value of prob
    indicates that a significant periodic signal is present.
  
  Reference: 
    Press, W. H. & Rybicki, G. B. 1989
    ApJ vol. 338, p. 277-280.
    Fast algorithm for spectral analysis of unevenly sampled data
    (1989ApJ...338..277P)

  Arguments:
      X   : Abscissas array, (e.g. an array of times).
      Y   : Ordinates array, (e.g. corresponding counts).
      Ofac : Oversampling factor.
      Hifac : Hifac * "average" Nyquist frequency = highest frequency
           for which values of the Lomb normalized periodogram will
           be calculated.
      
   Returns:
      Wk1 : An array of Lomb periodogram frequencies.
      Wk2 : An array of corresponding values of the Lomb periodogram.
      Nout : Wk1 & Wk2 dimensions (number of calculated frequencies)
      Jmax : The array index corresponding to the MAX( Wk2 ).
      Prob : False Alarm Probability of the largest Periodogram value
      MACC : Number of interpolation points per 1/4 cycle
            of highest frequency

  History:
    02/23/2009, v1.0, MF
      Translation of IDL code (orig. Numerical recipies)
  """
  #Check dimensions of input arrays
  n = long(len(x))
  if n != len(y):
    print 'Incompatible arrays.'
    return

  nout  = 0.5*ofac*hifac*n
  nfreqt = long(ofac*hifac*n*MACC)   #Size the FFT as next power
  nfreq = 64L             # of 2 above nfreqt.

  while nfreq < nfreqt: 
    nfreq = 2*nfreq

  ndim = long(2*nfreq)
  
  #Compute the mean, variance
  ave = y.mean()
  ##sample variance because the divisor is N-1
  var = ((y-y.mean())**2).sum()/(len(y)-1) 
  # and range of the data.
  xmin = x.min()
  xmax = x.max()
  xdif = xmax-xmin

  #extirpolate the data into the workspaces
  wk1 = zeros(ndim, dtype='complex')
  wk2 = zeros(ndim, dtype='complex')

  fac  = ndim/(xdif*ofac)
  fndim = ndim
  ck  = ((x-xmin)*fac) % fndim
  ckk  = (2.0*ck) % fndim

  for j in range(0L, n):
    __spread__(y[j]-ave,wk1,ndim,ck[j],MACC)
    __spread__(1.0,wk2,ndim,ckk[j],MACC)

  #Take the Fast Fourier Transforms
  wk1 = ifft( wk1 )*len(wk1)
  wk2 = ifft( wk2 )*len(wk1)

  wk1 = wk1[1:nout+1]
  wk2 = wk2[1:nout+1]
  rwk1 = wk1.real
  iwk1 = wk1.imag
  rwk2 = wk2.real
  iwk2 = wk2.imag
  
  df  = 1.0/(xdif*ofac)
  
  #Compute the Lomb value for each frequency
  hypo2 = 2.0 * abs( wk2 )
  hc2wt = rwk2/hypo2
  hs2wt = iwk2/hypo2

  cwt  = sqrt(0.5+hc2wt)
  swt  = sign(hs2wt)*(sqrt(0.5-hc2wt))
  den  = 0.5*n+hc2wt*rwk2+hs2wt*iwk2
  cterm = (cwt*rwk1+swt*iwk1)**2./den
  sterm = (cwt*iwk1-swt*rwk1)**2./(n-den)

  wk1 = df*(arange(nout, dtype='float')+1.)
  wk2 = (cterm+sterm)/(2.0*var)
  pmax = wk2.max()
  jmax = wk2.argmax()


  #Significance estimation
  #expy = exp(-wk2)          
  #effm = 2.0*(nout)/ofac       
  #sig = effm*expy
  #ind = (sig > 0.01).nonzero()
  #sig[ind] = 1.0-(1.0-expy[ind])**effm

  #Estimate significance of largest peak value
  expy = exp(-pmax)          
  effm = 2.0*(nout)/ofac       
  prob = effm*expy

  if prob > 0.01: 
    prob = 1.0-(1.0-expy)**effm

  return wk1,wk2,nout,jmax,prob


def getSignificance(wk1, wk2, nout, ofac):
  """ returns the peak false alarm probabilities
  Hence the lower is the probability and the more significant is the peak
  """
  expy = exp(-wk2)          
  effm = 2.0*(nout)/ofac       
  sig = effm*expy
  ind = (sig > 0.01).nonzero()
  sig[ind] = 1.0-(1.0-expy[ind])**effm
  return sig


def ReadStar(file,err_yn):
	
	if err_yn == 1:
		time,flux,err,sky=np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
		return time,flux,err,sky
	
	if err_yn != 1:
		time,flux=np.loadtxt(file,usecols=[0,1],unpack=True)
		return time,flux
		
		
def MakeFitsTable(x,y):
	
	from pyfits import Column
	
	c1=Column(name='x', format='E', array=x)
	c2=Column(name='y', format='E', array=y)
	
	tbhdu=pf.new_table([c1,c2])
	#print tbhdu.header.ascardlist() 

	name ='StarPosXY.fits'
	
	if os.path.isfile(name) == True:
		os.system('rm -rf StarPosXY.fits')
		
	tbhdu.writeto(name)
	
	return 0


def WCS_xy2rd():
	
	if os.path.isfile('StarPosRaDec.fits') == True:
		os.system('rm -rf StarPosRaDec.fits')
	
	os.system('wcs-xy2rd -w /Volumes/DATA/nites/Results/M71/M71Solved/tpl.wcs -i StarPosXY.fits -o StarPosRaDec.fits')

	return 0
	

def GetRaDec():
	
	RA,DEC=[],[]
	
	t=pf.open('StarPosRaDec.fits')
	tbdata=t[1].data
	
	for i in range(0,len(tbdata)):
		ra=tbdata[i][0]
		dec=tbdata[i][1]
	
		ra1=(ra/15)
		ra2=(fmod(ra1,1)*60)
		ra3=(fmod(ra2,1)*60)
		
		if len(str(ra3).split('.')[0]) < 2:
			ra3="0"+str(ra3)
		
		ratot="%02d:%02d:%s" % (int(ra1),int(ra2),str(ra3)[:5])
		RA.append(ratot)
		
		
		dec1=dec
		dec2=(fmod(dec1,1)*60)
		dec3=(fmod(dec2,1)*60)

		if len(str(dec3).split('.')[0]) < 2:
			dec3="0"+str(dec3)
		
		dectot="%02d:%02d:%s" % (dec1,dec2,str(dec3)[:5])
		DEC.append(dectot)
		
		
	return RA,DEC


def GetPos(starid,section,p):
	
	command='/Volumes/DATA/nites/Results/M71/Photometry/Batch1/A8.0D12S8C12/%s/M71_%s.coo' % (section, section)
	
	f=open(command).readlines()
	
	for i in range(0,len(f)):
		num,x,y=f[i].split()
		if int(num) == starid:
			
			xk=np.empty(1)
			yk=np.empty(1)
			
			# keep x and y
			xk[0]=float(x)
			yk[0]=float(y)

			imx=x
			imy=y
			
			
			if p > 0:
				print "X Y: %.2f %.2f" % (xk[0],yk[0])
	
	# get image subsection coords into tpl full 
	# frame coords for RA and DEC calculations
	if section == '1_1':
		# x(tplr1_1) = x(tpl) - 20 ; y(tplr1_1) = y(tpl) - 20
		xk[0]=xk[0]-20.0
		yk[0]=yk[0]-20.0
		if p > 0:
			print "x(tp1): %.2f y(tpl): %.2f" % (xk[0],yk[0])
	
	if section == '1_2':
		# x(tplr1_2) = x(tpl) - 20 ; y(tplr1_2) = y(tpl) + 492
		xk[0]=xk[0]-20.0
		yk[0]=yk[0]+492.0
		if p > 0:
			print "x(tp1): %.2f y(tpl): %.2f" % (xk[0],yk[0])

	if section == '2_1':
		# x(tplr1_2) = x(tpl) + 492 ; y(tplr1_2) = y(tpl) - 20
		xk[0]=xk[0]+492.0
		yk[0]=yk[0]-20.0
		if p > 0:
			print "x(tp1): %.2f y(tpl): %.2f" % (xk[0],yk[0])
	
	if section == '2_2':
		# x(tplr1_2) = x(tpl) + 492 ; y(tplr1_2) = y(tpl) + 492
		xk[0]=xk[0]+492.0
		yk[0]=yk[0]+492.0
		if p > 0:
			print "x(tp1): %.2f y(tpl): %.2f" % (xk[0],yk[0])	
	
	d1=MakeFitsTable(xk,yk)
	if d1 != 0:
		print "Problem making .fits table, exiting!"
		sys.exit()
	
	d2=WCS_xy2rd()
	if d2 != 0:
		print "Problem making .fits table, exiting!"
		sys.exit()
	
	ra,dec=GetRaDec()
	
	for i in range(0,len(ra)):	
		if p > 0:
			print "\n%s+%s" % (ra[i],dec[i])
	
	return ra,dec,imx,imy


def GetFlux(section):
	
	command='/Volumes/DATA/nites/Results/M71/Photometry/Batch1/A8.0D12S8C12/%s/M71_%s.flux' % (section, section)
	
	tflux,tfluxerr=np.loadtxt(command,usecols=[3,4],unpack=True)

	return tflux,tfluxerr


def GetVariablesAndRMS(templist,tflux):
	
	stddev=np.empty(len(templist))
	rms=np.empty(len(templist))
	plot_flux=np.empty(len(templist))
	
	# 1 night
	stddev_1=np.empty(len(templist))
	rms_1=np.empty(len(templist))

	nframes=[]
	
	# get the stddev for all
	for i in range(0,len(templist)):
		time,flux,err,sky=ReadStar(templist[i])

		# add the template flux from Sextractor	
		starid=int(templist[i].split('_')[1])
		flux=flux+tflux[(starid-1)]
		stddev[i]=std(flux)
		
		# get stddev for only data on 2011-08-04		
		c1=np.empty(len(time))
		c2=np.empty(len(time))
		for j in range(0,len(time)):
			c1[j]=abs(time[j]-2455778.39028)
			c2[j]=abs(time[j]-2455778.67986)
		
		start=np.where(c1==min(c1))[0][0]
		end=np.where(c2==min(c2))[0][0]
				
		# get the fractional rms for all nights
		# get the plot flux 
		rms[i]=stddev[i]/tflux[(starid-1)]
		plot_flux[i]=tflux[(starid-1)]
		
		# get the fractional rms for 2011-08-04
		stddev_1[i]=std(flux[start:end])
		rms_1[i]=stddev_1[i]/tflux[(starid-1)]
		
		# only run this part when needed, typically once per batch
		get_nightly=0
		if get_nightly > 0:
			# try getting the nightly stddev per lc
			rms_nightly=GetNightlyRMS(templist,time,flux)
		if get_nightly==0:
			rms_nightly=0
		
		nframes.append(len(flux))
		print "[GetVariables] %d/%d" % ((i+1),len(templist))
		
	# variables	for investigation
	n=np.where(stddev>(2*(median(stddev))))
	# non-variables for RMS vs Flux plots
	n2=np.where(stddev<(2*(median(stddev))))
	
	# limit line for plot
	limitx=[0,(len(templist))]
	limity=[(median(stddev)*2),(median(stddev)*2)]
	
	figure(2)
	plot(stddev,'k-')
	plot(limitx,limity,'r-')
	ylabel('stddev')
	xlabel('image number')
	ylim(0,5000)
	
	
	figure(3)
	semilogx(plot_flux[n2],rms_1[n2],'r.')
	ylabel('rms (mag)')
	xlabel('flux')
	
	show()
	
	return n,n2,nframes,rms,rms_1,rms_nightly,plot_flux


def GetNightlyRMS(templist,time,flux):

	rms_nightly=np.empty((len(templist),43))

	start=[0]
	end=[]
	for j in range(0,len(time)-1):
		if abs(time[j+1]-time[j]) > 0.5:
			end.append(j)
			start.append(j+1)
	
	end.append(len(time)-1)		
	
	for j in range(0,len(start)):
		rms_nightly[i,j]=(std(flux[start[j]:end[j]]))/tflux[(starid-1)]	

	return rms_nightly


def GetBestNight(rms_nightly):
	loc=[]
	for i in range(0,len(rms_nightly)):
		loc.append(np.where(rms_nightly[i]==min(rms_nightly[i]))[0][0])

	d={}
	for i in set(loc):
		d[i]=loc.count(i)

	return d
	
	
def GetLombPeaks(fy,fx):

	lombp=zip(fy,fx)
	lombp.sort()
	fysort,fxsort=zip(*lombp)
	
	# makes a list of values, but they are backwards
	# i.e. best period is last
	peaks=list(fysort[-20:])
	freqs=list(fxsort[-20:])
	
	# reverse the values before returning them
	peaks.reverse()
	freqs.reverse()
		
	return peaks, freqs


def MakeSimbadList(n,templist,section):

	# name the file
	name="simlist_%s.txt" % (section)

	# check if it exists
	# remove it if it does
	if os.path.isfile(name) == True:
		print "\nOverwritting old simlist file..."
		command='rm -rf %s' % (name)
		os.system(command)
	
	# write out the data to file
	f=open(name,'w')
	
	for i in range(0,len(n[0])):
		vstar=templist[n[0][i]]
		starid=int(vstar.split('_')[1])		

		ra,dec,imx,imy=GetPos(starid,section,0)
		line="%s+%s\n" % (ra[0],dec[0])
		f.write(line)
		
	f.close()	
	
	return 0


def MakeRegions(section,imx,imy,starid):

	# Add new star to region file
	# template region : circle(36.88,481.42,8) # text={94}
	file='/Users/James/Documents/Observing/NITESObs/M71/Batch1/%s/%sVariables.reg' % (section,section)
	
	# set up the file is not there already
	if os.path.isfile(file) == False:
		command='touch /Users/James/Documents/Observing/NITESObs/M71/Batch1/%s/%sVariables.reg' % (section,section)
		os.system(command)
		mf=open(file,'a')
		templine1='# Region file format: DS9 version 4.1\n'
		templine2='# Filename: tplr%s.fits\n' % (section)
		templine3='global color=green dashlist=8 3 width=1 font="helvetica 10 normal" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
		templine4='physical\n'
		
		mf.write(templine1)
		mf.write(templine2)
		mf.write(templine3)
		mf.write(templine4)
	
		mf.close()
	
	f=open(file,'a')
	
	line="circle(%.2f,%.2f,8) # text={%d}\n" % (float(imx),float(imy),starid)
	f.write(line)
	f.close()
	
	return 0


def GetBinnedLc(time_n,flux,err,sxphe):
	
	# zip lists and sort them in order of phase
	# unzip for binning
	temp=zip(time_n,flux,err)
	temp.sort()
	tsort,fsort,errsort=zip(*temp)
	
	if sxphe < 1.0:
		bin=1000.0
	
	if sxphe > 0.0:
		bin=len(flux)/3
	
	# binning factor to give 'bin' points in final lc
	bf=int(len(temp)/bin)
	
	# binned arrays
	tbinned=np.empty(bin)
	fbinned=np.empty(bin)
	errbinned=np.empty(bin)
	
	for j in range(0,len(tbinned)):
		tbinned[j]=average(tsort[(j*bf):bf*(j+1)])
		fbinned[j]=average(fsort[(j*bf):bf*(j+1)])
		errbinned[j]=average(errsort[(j*bf):bf*(j+1)])
		
	return tbinned, fbinned, errbinned


def GetMags(flux):

	mags=np.empty(len(flux))

	for i in range(0,len(flux)):
		mags[i]=(-2.5*log10(flux[i])) + 25.0

	return mags


def GetTransit(tbinned,fbinned):

	pvals=np.empty(len(tbinned))
	transit_new=np.empty(len(tbinned))
	transit_mag=np.empty(len(tbinned))

	# masked array for low order polynomial fit   
	index=np.where((tbinned<0.2)+(tbinned>0.8)) 
	
	print "1st or 2nd Order Fit?"
	fittype = raw_input("(1/2): ")
	
	if str(fittype) == '1':
	
		# best line of fit for flux out of transit
		coeffs=polyfit(tbinned[index],fbinned[index],1)
	
		# best vals for mags
		besty=polyval(coeffs,tbinned[index])
	
		# finding p to fit the data -- equation -- p2xk^2 + p1xk + p0 =yk
		for i in range(0,len(fbinned)):
			pvals[i]=coeffs[1]+(coeffs[0]*tbinned[i])

		for i in range(0,len(fbinned)):
			transit_new[i]=fbinned[i]/pvals[i]
	
	
	if str(fittype) == '2':
	
		# best line of fit for flux out of transit
		coeffs=polyfit(tbinned[index],fbinned[index],2)
	
		# best vals for mags
		besty=polyval(coeffs,tbinned[index])
	
		# finding p to fit the data -- equation -- p2xk^2 + p1xk + p0 =yk
		for i in range(0,len(fbinned)):
			pvals[i]=coeffs[2]+(coeffs[1]*tbinned[i])+(coeffs[0]*(tbinned[i]*tbinned[i]))

		for i in range(0,len(fbinned)):
			transit_new[i]=fbinned[i]/pvals[i]
	
	
	# remember to do errors later
	for i in range(0,len(fbinned)):
		transit_mag[i]=-2.5*log10(transit_new[i])
	
	
	figure(20)	
	plot(tbinned,fbinned,'r.')
	plot(tbinned[index],besty,'-k')
	
	figure(21)
	plot(tbinned,transit_new,'r.')
	
	show()
	
	return transit_mag


# print phased light curve for plotting in gnuplot 
def PrepareTamuz(templist,nframes):
	
	if os.path.exists('tamuz') == False:
		os.mkdir('tamuz')
		os.chdir('tamuz')
	if os.path.exists('tamuz') == True:
		os.chdir('tamuz')
	
	for i in range(0,len(templist)):
		if nframes[i]==25741:
			if os.path.isfile(templist[i]) == False:
				comm="cp ../%s ." % (templist[i])
				os.system(comm)
				
	return 0


def PrintFile(section,starid,time,time2,flux,err,period,bub):
	
	savedir="/Volumes/DATA/nites/Results/M71/Photometry/Batch1/Variables/%s/final" % (section)
	
	if os.path.exists(savedir) == False:
		os.mkdir(savedir)
	
	if bub == 1:
		name = "%s/star_%05d_%s_b_FIN.lc.txt" % (savedir,starid, section)
	if bub == 0:
		name = "%s/star_%05d_%s_ub_FIN.lc.txt" % (savedir,starid, section)
	
	t_out=np.concatenate((time,time2))
	f_out=np.concatenate((flux,flux))
	err_out=np.concatenate((err,err))
	
	z=np.concatenate((t_out,f_out,err_out)).reshape(3,len(t_out)).transpose()
	np.savetxt(name,z,fmt='%.8f    %.5f    %.5f')
	
	# open the file and get contents
	f=open(name,'r')
	s=f.readlines()
	f.close()
	
	# reopen in to write title line and contents back
	f=open(name,'w')
	line="# Star: %05d Period: %.6f\n" % (starid,period)
	f.write(line)
	for i in range(0,len(s)):
		f.write(s[i])
	f.close()
	
	return 0


def Normalize(starid,time,time2,flux,err,period,bub):
	
	nf=np.copy(flux)
	nf.sort()
	
	norm_f=average(nf[-100:])
	
	f_norm=flux/norm_f
	f_mag=-2.5*log10(f_norm)
	f_magerr=err/flux[0]
	
	figure(10)
	errorbar(time,f_mag,yerr=f_magerr,fmt='r.')
	errorbar(time2,f_mag,yerr=f_magerr,fmt='r.')
	gca().invert_yaxis()
	title("Star %d (binned)" % (starid))
	ylabel("Differential Magnitude")
	xlabel("Phase")
	
	show()
	
	# work out V_max and dV after normalising
	# 26.974 was worked out using AH 1971 standards
	# see notes
	Vmax=(-2.5*log10(norm_f))+26.974

	dV=GetDeltaV(f_mag,bub)

	print "V_max = %.2f [%.2f]" % (Vmax, norm_f)
	print "dV = %.4f" % (dV)
	
	print "Final Output?"
	yn=raw_input("(e.g. y): ")
	if str(yn) == 'y':
		output=PrintFile(section,starid,time,time2,f_mag,f_magerr,period,bub)
		if output != 0:
			print "Problem printing output file, exiting!"
			sys.exit()

	return f_norm,f_mag,f_magerr


def GetPeriodError(section,starid):

	file="/Volumes/DATA/nites/Results/M71/Photometry/Batch1/A8.0D12S8C12/%s/period/star_%05d_%s_logfile.dat" % (section,starid,section)
	
	f=open(file).readlines()
	
	period_err=float(f[-2].split()[-1])
	
	return period_err


def GetDeltaV(f_mag,bub):
	
	fmag_cp=np.copy(f_mag)
	
	fmag_cp.sort()
	
	if bub == 0:
		v_min = average(fmag_cp[-100:])
		v_max = average(fmag_cp[:100])
	
	if bub == 1:
		v_min = average(fmag_cp[-10:])
		v_max = average(fmag_cp[:10])
	
	dv=v_min-v_max
	
	return dv


###############################################
################### MAIN ######################
###############################################

# toggles
var_rms = 0
make_simlist=0
prep_tamuz=0
run_p = 1
sxphe = 0
err_yn = 1


# get image list
templist=commands.getoutput('ls star_0016*.lc.txt').split('\n')
# get image subsection
section="%s_%s" % (templist[0].split('_')[2],templist[0].split('_')[3][0])

# get the sextractor fluxes
tflux,tfluxerr=GetFlux(section)

# find the variables and RMS's etc
if var_rms > 0:
	n,n2,nframes,rms,rms_1,rms_nightly,plot_flux=GetVariablesAndRMS(templist,tflux)

	if rms_nightly != 0:
		d=GetBestNight(rms_nightly[n2])

n=(array([0]),)

# or use previous variables 
#n=(array([0, 4, 5, 6, 47, 48, 49, 50, 54, 57, 58, 59, 84, 90, 91, 94, 119, 122, 123, 181, 182, 194, 206, 223, 228, 250, 251, 277, 280, 282, 283, 304, 321, 330, 362, 365, 377, 380, 398, 420, 436, 445, 447, 453, 467, 472, 475, 486, 491, 493, 494, 497, 502, 506, 508, 518, 522, 533, 548, 556, 604, 630, 633, 634, 635, 642, 662, 664, 665, 670, 677, 704, 714, 715, 717, 818, 824, 843, 851, 855, 896, 924, 930, 937]),) # M71 b1 1_1 

# n=(array([2, 3, 5, 14, 16, 17, 29, 56, 81, 83, 102, 105, 121, 129, 130, 131, 133, 136, 139, 141, 150, 151, 157, 159, 174, 175, 185, 248, 251, 261, 262, 270, 274, 281, 285, 286, 287, 288, 290, 292, 296, 304, 305, 306, 307, 313, 314, 316, 323, 327, 336, 338, 342, 356, 363, 365, 368, 369, 390, 398, 411, 416, 417, 420, 423, 440, 445, 446, 457, 462, 463, 467, 472, 473, 474, 485, 494, 495, 502, 508, 512, 535, 538, 565, 566, 573, 574, 583, 600, 618, 637, 676, 688, 689, 693, 694, 705, 737, 751, 770, 771, 773, 775, 776, 777, 778, 779, 803, 806, 849, 865, 928, 935, 936, 1012, 1028, 1036, 1037, 1041, 1046, 1066, 1068, 1071, 1073, 1076, 1077, 1081, 1086, 1087, 1089, 1090, 1103, 1104]),) # M71 b1 1_2

# n=(array([0, 3, 4, 5, 6, 7, 8, 9, 10, 11, 14, 18, 21, 23, 24, 33, 37, 39, 40, 43, 44, 52, 62, 72, 79, 82, 87, 98, 99, 100, 108, 117, 118, 122, 130, 131, 138, 139, 148, 153, 159, 165, 167, 170, 171, 179, 180, 187, 196, 197, 198, 199, 208, 209, 210, 211, 214, 223, 236, 241, 242, 243, 244, 251, 314, 315, 322, 327, 331, 346, 348, 350, 351, 355, 367, 368, 369, 370, 373, 381, 383, 388, 389, 390, 391, 396, 402, 409, 415, 420, 421, 427, 445, 473, 488, 490, 492, 520, 527, 579, 591, 597, 603, 612, 614, 618, 621, 622, 623, 626, 632, 638, 666, 674, 678, 695, 698, 711, 717, 725, 726, 728, 743, 745, 748, 752, 753, 754, 762, 764, 765, 781, 798, 799, 807, 815, 849, 863, 872, 878, 892, 902, 906, 910, 914, 922, 931, 959]),) # M71 b1 2_1

# n=(array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 13, 15, 27, 30, 63, 185, 218, 232, 233, 245, 273, 280, 298, 300, 319, 323, 328, 329, 333, 334, 342, 346, 355, 356, 359, 360, 361, 362, 363, 364, 365, 369, 372, 376, 382, 386, 387, 391, 396, 397, 399, 402, 406, 407, 408, 417, 418, 426, 430, 436, 442, 443, 452, 455, 456, 457, 458, 459, 463, 467, 472, 482, 486, 498, 502, 503, 523, 537, 541, 578, 579, 648, 694, 702, 716, 724, 759, 761, 769, 771, 772, 775, 776, 777, 778, 786, 792, 856, 871, 889, 893, 894, 911, 912, 917, 921, 922, 923, 931, 941, 953, 986]),) # M71 b1 2_2

# Make a list for simbad?
if make_simlist > 0:
	simlist=MakeSimbadList(n,templist,section)
	if simlist != 0:
		print "Problem making Simbad list, exiting..."
		sys.exit()

# Prepare files for Tamuz
if prep_tamuz > 0:
	done=PrepareTamuz(templist,nframes)
	if done != 0:
		print "Problem making preparing for Tamuz, exiting..."
		sys.exit()

# Run period + phasing loop
if run_p > 0:
	
	# now run lomb-scargle on the variables
	for i in range(0,len(n[0])):
		
		vstar=templist[n[0][i]]
		
		if err_yn == 1:
			time,flux,err,sky=ReadStar(vstar,err_yn)
		
		if err_yn != 1:
			time,flux=ReadStar(vstar,err_yn)
			
		starid=int(vstar.split('_')[1])
	
		# add the template flux from Sextractor
		if sxphe < 1.0:
			flux=flux+tflux[(starid-1)]
		
		# run lomb-scargle
		fx,fy, nout, jmax, prob = fasper(time,flux, 6., 6.)
		
		# get error from scargle
		period_err=GetPeriodError(section,starid)
		
		# get top 20 peaks and periods
		peaks,freqs=GetLombPeaks(fy,fx)
		
		# get this as a check
		ls_period=1/fx[jmax]
		
		figure(4)
		title('Lomb-Scargle')
		ylabel('Power')
		xlabel('Period')
		plot((1/fx),fy,'r-')
		xlim(0,50)
		xticks(arange(0,52,2))
		
		
		print "\n----------------------------------------------------------"
		print "Star:\t\t\t%s [%d/%d]" % (templist[n[0][i]],(i+1),len(n[0]))
		print "Template Flux:\t\t%f" % (tflux[(starid-1)])
		print "Template Err:\t\t%f" % (tfluxerr[(starid-1)])
		print "Number of Points:\t%d" % (len(flux))	
		for j in range (0,len(peaks)):
			print "\tPeriod[%d]:\t%f" % ((j+1),(1/(freqs[j])))
		#print "Removed %d exteme outliers" % (n_outliers)
		print "----------------------------------------------------------\n"
		print "LS Period: %f (%f)" % (ls_period,period_err)
		print "Check periods!"
	
		choice = 0.0
		run=0
		redo=0
		while choice < 1.0:
			
			# redo = 1 if only running 'gs'
			# no need to plot everything again when only getting the coordinates
			if redo < 1.0:	
				epochs=zeros(2000,float)	
				
				epoch=time[0]
				
				if run < 1.0:
					period=ls_period
					#ra,dec,imx,imy=GetPos(starid,section,1)
					# making regions is done!
					#marked=MakeRegions(section,imx,imy,starid)
					print "Star ID: %d\n" % (starid)
				
				if run > 0.0:
					period=pnew
					print "New Period: %f" % (period)
				
				time_n=np.empty(len(time))
				
				for j in range(0,len(time)):
					
					x=((time[j]-time[0])/period)%1.0
					time_n[j]=x
				
				
				if err_yn == 1:
					# get binned lc
					tbinned,fbinned,errbinned=GetBinnedLc(time_n,flux,err,sxphe)
				
					# set min light to phase 0
					loc=np.where(fbinned==min(fbinned))
					tbinned=tbinned-tbinned[loc]
				
					# double the number of cycles for plotting
					tbinned2=np.empty(len(tbinned))
					time_n2=np.empty(len(time_n))
				
					for j in range(0,len(tbinned)):
						tbinned2[j]=tbinned[j]+1.0
				
					for j in range(0,len(time_n)):
						time_n2[j]=time_n[j]+1.0	
				
				
					# get mags for plots
					fluxm=GetMags(flux)
					fbinnedm=GetMags(fbinned)
				
				figure(5)	
				plot(time_n,flux,'r.')
				if err_yn == 1:
					plot(time_n2,flux,'r.')
				title('Star %d: Raw Data (P=%.5f)' % (starid,period))
				xlabel('Phase')
				ylabel('Flux')
				xlim(0,2.0)
				
				if err_yn == 1:
					figure(6)
					plot(tbinned,fbinned,'r.')
					plot(tbinned2,fbinned,'r.')
					title('Star %d: Binned Data (P=%.5f)' % (starid,period))
					xlabel('Phase')
					ylabel('Flux')
					#xlim(0,2.0)
				show()
			
			
			line=raw_input("")
			
			if str(line[0]) == '0':
				choice=float(line.split()[0])
				pnew=float(line.split()[1])
				run = run + 1
	
			if str(line[0]) == '1':
				choice=1.0
				break
			
			if str(line) == 'ft':
				transit_mag=GetTransit(tbinned,fbinned)
				redo = 1.0
				choice = 0.0
			
			if str(line) == 'pfb':
				output=PrintFile(section,starid,tbinned,tbinned2,fbinned,errbinned,period,1)
				if output != 0:
					print "Problem printing output file, exiting!"
					sys.exit()
				
				choice = 1.0
			
			if str(line) == 'pfub':
				output=PrintFile(section,starid,time_n,time_n2,flux,err,period,0)
				if output != 0:
					print "Problem printing output file, exiting!"
					sys.exit()
				
				choice = 1.0
			
			if str(line) == 'normb':
				f_norm,f_mag,f_magerr=Normalize(starid,tbinned,tbinned2,fbinned,errbinned,period,1)
				choice = 1.0
			
			if str(line) == 'normub':
				f_norm,f_mag,f_magerr=Normalize(starid,time_n,time_n2,flux,err,period,0)
				choice = 1.0