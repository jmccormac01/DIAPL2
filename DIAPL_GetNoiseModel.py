
# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_GetNoiseModel.py - 	a program to get the Noise Model after DIAPL2 reduction
#						
# 						
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 23/11/11 - 	code writen
# 24/11/11 -	code tested
# 28/11/11 - 	added rms vs bin code from Tamuz_New.py 
#


import commands as cmd
import numpy as np
import pyfits as pf
import os, sys
import time
import pylab as pl
from pyraf import iraf

def GetSkyMedian(gain,trim):
	
	# for this and subtracted images there is a new way to get the sky level.
	# this involves measuring the stddev in the subtracted images and then squaring it 
	# and multiplying by the gain to get the flux per pix which can then be used in the
	# sky calculations, see below to do so.
	
	# first imcopy the subtracted images into a new file
	# this is to stopp the truncation problem which is common. 
	
	t=cmd.getoutput('ls /Volumes/DATA/NITES/data/2011-08-12/reduced/CheckSkyInSubtracted/s_M71*.fits').split('\n')
	
	if trim == 1:
		for i in range(0,len(t)):
			image=str(t[i]+"[1:500,1:500]")
			image2=str(t[i])

			# clobber the old image
			iraf.imcopy(input=image,output=image2)
	
	# then define a box in each image (same box for bright and dark time) and 
	# get the standard deviation in it for each frame that night. Then get the average 
	# nightly stddev which is used to get the sky level.
	
	std_list=np.empty(len(t))
	
	for i in range(0,len(t)):
		h=pf.open(t[i])
		d=h[0].data[243:313,274:344]
		
		std_list[i]=np.std(d)
		#print "%s std: %.2f" % (t[i],std_list[i])
		
	av=pow(np.average(std_list),2)*gain
	print "Sky_bg_pp: %.2f e-" % (av)
	
	
	return av
	
	
def FindCalibs():
	
	odir=os.getcwd()
	
	print "\nFinding calibs..."
	stop=10
	while stop > 1:
		time.sleep(2)
		pwd=os.getcwd()
		print "Searching " + str(pwd)
		templist=os.listdir(pwd)
		
		for i in range(0,len(templist)):
			if str(templist[i]) == 'calibs':
				path=str(pwd)+'/calibs/'
				stop=0
				print "Calibs found at " + str(path)	
				print "Returning to " + str(odir)
				os.chdir(odir)
				return path
		os.chdir('../')
		
	return 0	


def GetTotalFlatCounts(path,gain):
	
	line1='ls %s/NITES-AutoFlat*' % (path)
	
	# get a list of the flat fields
	templist=cmd.getoutput(line1).split('\n')

	fmean=[]
	
	for i in range(0,len(templist)):
		hdulist=pf.open(templist[i])
		data=hdulist[0].data
		
		fmean.append(np.mean(data)*gain)
		
	total_flux=np.sum(fmean)
	err=np.sqrt(total_flux)
	
	return total_flux, err


def GetDarkCurrent(path,gain):
	
	line=str(path) + 'Dark.fits'
	
	hdulist=pf.open(line)
	data=hdulist[0].data

	exptime=float(hdulist[0].header['EXPTIME'])
		
	fmean=np.mean(data)*gain
		
	dc=fmean/exptime
	
	return dc


def ReadStar(file):
	
	s=open(file,'r').readlines()
	
	time=np.empty(len(s))
	flux=np.empty(len(s))
	
	for i in range(0,len(s)):
		time[i]=float(s[i].split('\n')[0].split()[0])
		flux[i]=float(s[i].split('\n')[0].split()[1])
		
	return time,flux


def GetFlux(section):
	
	command='/Volumes/DATA/NITES/Results/M71/Photometry/Batch1/A8.0D12S8C12/%s/M71_%s.flux' % (section, section)
	
	f=open(command).readlines()
	
	tflux=np.empty(len(f))
	tfluxerr=np.empty(len(f))
	
	for i in range(0,len(f)):
		tflux[i]=float(f[i].split('\n')[0].split()[3])
		tfluxerr[i]=float(f[i].split('\n')[0].split()[4])

	return tflux,tfluxerr
	

def GetRMSvsFlux(section,tflux):
	
	comm='ls /Volumes/DATA/NITES/data/2011-08-12/reduced/phot/A6/%s/star_0*.lc.txt' % (section)
	#comm='ls /Volumes/DATA/NITES/Results/M71/Photometry/Batch1/BLS/%s/binned/2011-08-04/star_0*.lc.txt' % (section)
	
	templist=cmd.getoutput(comm).split('\n')

	stddev=np.empty(len(templist))
	rms=np.empty(len(templist))
	plot_flux=np.empty(len(templist))
	
	starid=np.empty(len(templist))
	
	# get the stddev for all
	for i in range(0,len(templist)):
		time,flux=ReadStar(templist[i])

		# add the template flux from Sextractor	
		starid[i]=int(templist[i].split('_')[2])
		stddev[i]=np.std(flux)
				
		plot_flux[i]=tflux[(starid[i]-1)]
		
		# get the fractional rms for 2011-08-04
		stddev[i]=np.std(flux)
		rms[i]=stddev[i]/tflux[(starid[i]-1)]
				
		#print "[RMS]: %d/%d" % ((i+1),len(templist))
		
	# non-variables for RMS vs Flux plots
	n=np.where(stddev<(2*(np.median(stddev))))
	n2=np.where(stddev>=(2*(np.median(stddev))))
	
	return n,n2,rms,plot_flux,starid
	
	
def GetNoiseModel(airmass,exptime,sky_bg_pp,gain): #,,data_fnm,bins): 
	
	# we need to get the error models in two formats:
	# 	1 - for rms vs magnitude plots
	# 	2 - for rms vs binsize plots
	# only those that change with exposure time will be included.
	# those are:
	#	1 - Object Noise
	#	2 - Sky Noise
	# 	3 - Dark Current
	#	4 - Scintillation
	
	
	#print "Phot aperture radius?"
	#Aper=float(raw_input("(e.g. 2.5): "))
	Aper=6.0
	
	#print "Airmass Exponent?"
	#print "1.50 Perpendicular to wind"
	#print "1.75 Close to Zenith"
	#print "2.00 Parallel to wind"
	#Aexp=float(raw_input("(e.g. 1.75): "))
	Aexp=1.75
	
	print "\nAssuming: "
	print "\tNITES ReadNoise: 10.0 e/pix"
	print "Calculating noise models..."
	
	# RMS VS MAGNITUDE PLOT ERRORS #	
	# list for model backgrounds
	sbg_model=[]
	object_model=[]
	rn_model=[]
	dc_model=[]
	ff_model=[]
	scint_model=[]
	total_err_model_wff=[]
	total_err_model_wscintff=[]
	
	# example star fluxes
	fluxes=np.linspace(1000,1000000,400)
	
	# read noise model
	Npix=np.pi*(Aper*Aper)
	#ReadNoise_pp=14.0
	ReadNoise_pp=10.0
	ReadNoise=ReadNoise_pp*Npix
	
	for i in range(0,len(fluxes)):
		rn_model.append(np.sqrt((ReadNoise_pp*ReadNoise_pp*Npix))/fluxes[i])
	
	# average error in sky background
	sky_bg=sky_bg_pp*Npix
	sbg_err=np.sqrt(sky_bg)
	# sky model
	for i in range(0,len(fluxes)):
		sbg_model.append(sbg_err/fluxes[i])
	
	
	path=FindCalibs()
	if path == 0:
		print "Could not find calibs, exiting!"
		sys.exit()
	
	# flat field model
	tf,err=GetTotalFlatCounts(path,gain)
	ff_frac_err=err/tf
		
	for i in range(0,len(fluxes)):
		ff_model.append(ff_frac_err)
	
	# dark current model
	dc_pp_ps=GetDarkCurrent(path,gain)
	dc_pp_ps=dc_pp_ps/2.0
	DarkCurrent=dc_pp_ps*exptime*Npix
	
	for i in range(0,len(fluxes)):
		dc_model.append(np.sqrt(DarkCurrent)/fluxes[i])
	
	# object model
	for i in range(0,len(fluxes)):
		object_model.append(np.sqrt(fluxes[i])/fluxes[i])	
	
	# summary so far
	print "\n"
	print "Total flux in Flat.fits:\t" + str(tf)
	print "Error in Flat.fits:\t\t" + str(err)
	print "Fractional error:\t\t" + str(ff_frac_err)
	print "Dark Current (ppps):\t\t" + str(dc_pp_ps)
	print "Total Dark Current (pape):\t" + str(DarkCurrent)
	
	# scintillation models
	# done for bin sizes too
	D=18.0
	ho=8000.0
	h=2332.0
	
	EXPTIME=np.linspace(30,12000,400)
	AvAirmass=airmass[len(airmass)/2]
	scint=np.empty(len(EXPTIME))

	w1=exptime*2
	
	# scintillation model for nominal exposure time
	for i in range(0,len(fluxes)):
		scint_model.append(0.09*(pow(D,(-2.0/3.0)))*(pow(AvAirmass,float(Aexp)))*(np.exp(-h/ho))*(pow(w1,(-1.0/2.0))))

	# total error from vik dhilons ccd equation + flats
	for i in range(0,len(fluxes)):
		total_err_model_wff.append(pow((fluxes[i]+sky_bg+DarkCurrent+((ReadNoise_pp*ReadNoise_pp)*Npix))+np.sqrt(fluxes[i]*Npix*ff_model[i]),0.5)/fluxes[i])
	
	#
	
	# add on scintillation
	for i in range(0,len(fluxes)):
		total_err_model_wscintff.append(np.sqrt((pow(total_err_model_wff[i],2))+(pow(scint_model[i],2))))
	
	return fluxes, sbg_model, object_model, rn_model, dc_model, ff_model, scint_model, total_err_model_wff, total_err_model_wscintff


def GetRMSvsBin(plot_flux,rms,starid,n,exptime,time_a,airmass):
	
	# analyze data based on brightness
	# only use the ones flagged by n as less variable
	temp=zip(plot_flux[n],rms[n],starid[n])
	
	# remove duplicate items
	# remnant of ngts data, might not apply here
	temp=list(set(temp))
	temp.sort()
	
	# unzip back to the arrays
	plot_flux_s,rms_s,starid_s=zip(*temp)
	
	# make them numpy arrays
	plot_flux_s=np.copy(plot_flux_s[::-1])
	rms_s=np.copy(rms_s[::-1])
	starid_s=np.copy(starid_s[::-1])
		
	
	# get top lim brightest star data
	if section == "1_1":
		lim=20
	if section == "1_2":
		lim=20
	if section == "2_1":
		lim=20
	if section == "2_2":
		lim=20

	for i in range(0,lim):
	
		# these flux values already have tflux added on!
		file="/Volumes/DATA/NITES/data/2011-08-04/reduced/phot/A6/%s/star_%05d_%s_2011-08-04.lc.txt" % (section,starid_s[i],section)
		
		if i == 0:
			time_t,flux_t=ReadStar(file)
		if i > 0:
			time,flux=ReadStar(file)
			time_t=np.vstack([time_t,time])
			flux_t=np.vstack([flux_t,flux])
	
	# define the bins for each star
	# bin until there are only 5 binned points
	z=1
	bins=[0]*(len(time_t[0])/5)
	for i in range(0,len(bins)):
		bins[i]=z
		z=z+1
		
	pvals=np.empty( (len(time_t),len(time_t[0])) )
	flux_new=np.empty( (len(time_t),len(time_t[0])) )

	# try normalising the curves with a polynomial
	# to solve the rms vs bin deviation
	for k in range(0,len(time_t)):
	
		# best line of fit for flux 
		coeffs=np.polyfit(time_t[k],flux_t[k],1)
		
		# best vals for mags
		besty=np.polyval(coeffs,time_t[k])
		
		# finding p to fit the data -- equation -- p2xk^2 + p1xk + p0 =yk
		for i in range(0,len(time_t[k])):
			#pvals[k][i]=coeffs[2]+(coeffs[1]*time_t[k][i])+(coeffs[0]*(time_t[k][i]*time_t[k][i]))
			pvals[k][i]=coeffs[1]+(coeffs[0]*time_t[k][i])
	
		for i in range(0,len(time_t[1])):
			flux_new[k][i]=(flux_t[k][i]/pvals[k][i])*np.average(flux_t[k])
			
		print "[Before] std: %f" % (np.std(flux_t[k])/np.average(flux_t[k]))
		print "[After] std: %f" % (np.std(flux_new[k])/np.average(flux_t[k]))
	
	
	# un-normalised data
	rms_d_b=np.empty((len(time_t),len(bins)))
	exp_time=np.empty((len(time_t),len(bins)))
			
	for k in range(0,len(time_t)):
		
		for j in range(0,len(bins)):
			d_b=np.empty(len(time_t[k])/bins[j])
		
			q=0
			for i in range(0,len(time_t[k])/bins[j]):
				d_b[i]=((np.sum(flux_t[k][q:q+bins[j]]))/bins[j])
				q=q+bins[j]
		
			# rms in mag
			rms_d_b[k][j]=(np.std(d_b)/np.average(flux_t[k]))
			# exp time in secs per bin
			exp_time[k][j]=bins[j]*exptime		
	
	
	# normalised data
	rms_d_b_n=np.empty((len(time_t),len(bins)))
	exp_time_n=np.empty((len(time_t),len(bins)))
			
	for k in range(0,len(time_t)):
		
		for j in range(0,len(bins)):
			d_b_n=np.empty(len(time_t[k])/bins[j])
		
			q=0
			for i in range(0,len(time_t[k])/bins[j]):
				d_b_n[i]=((np.sum(flux_new[k][q:q+bins[j]]))/bins[j])
				q=q+bins[j]
		
			# rms in mag
			rms_d_b_n[k][j]=(np.std(d_b_n)/np.average(flux_t[k]))
			# exp time in secs per bin
			exp_time_n[k][j]=bins[j]*exptime	
	
	
	# now median out the rms vs bin size from the 9 stars
	rms_med=np.copy(rms_d_b.transpose())
	rms_med_n=np.copy(rms_d_b_n.transpose())
	#rms_med_final=np.median(rms_med,axis=1)
	
	# stars seemed to be varying, only 1 + 7 showed nice plots - use them only
	rms_med_final=(rms_med[:,0]+rms_med[:,6])/2
	# AFTER NORMALISATION 
	# 1_1 batch 1 - stars 1,2,4,5 + 7 showed nice plots - use them only
	if section == "1_1":
		rms_med_final_n=(rms_med_n[:,0]+rms_med_n[:,1]+rms_med_n[:,3]+rms_med_n[:,4]+rms_med_n[:,6])/5
	if section == "1_2":
		rms_med_final_n=(rms_med_n[:,2]+rms_med_n[:,5]+rms_med_n[:,8]+rms_med_n[:,12])/4
	if section == "2_1":
		rms_med_final_n=(rms_med_n[:,6]+rms_med_n[:,13]+rms_med_n[:,15])/3
	if section == "2_2":
		rms_med_final_n=(rms_med_n[:,1]+rms_med_n[:,2]+rms_med_n[:,4]+rms_med_n[:,6]+rms_med_n[:,9]+rms_med_n[:,11])/6
	
	# now workout the 1/N^0.5 line for plotting gausian noise
	bins_longer=np.arange(1,201,1)
	binsr=(1/np.sqrt(bins_longer))*rms_med_final[0]
	
	binsr_n=(1/np.sqrt(bins_longer))*rms_med_final_n[0]
	
	# look at tms vs bin plots for top stars?
	cycle_top = 1
	if cycle_top > 0:
	
		for j in range(0,len(rms_d_b)):
			bins_longer=np.arange(1,201,1)
			binsr=(1/np.sqrt(bins_longer))*rms_d_b[j][0]	
			
			binsr_n=(1/np.sqrt(bins_longer))*rms_d_b_n[j][0]
			
			
			pl.figure(j+2)
			pl.title("star %d [%d]" % ((j+1),starid_s[j]))
			pl.loglog(exp_time[0],rms_d_b[j],'r.',(bins_longer)*exptime,binsr,'k--')
			pl.loglog(exp_time[0],rms_d_b_n[j],'g.',(bins_longer)*exptime,binsr_n,'b--')
			
			pl.figure(20+j)
			pl.subplot(211)
			pl.title("star %d [%d]" % ((j+1),starid_s[j]))
			pl.plot(time_t[j],flux_t[j],'r.')
			pl.subplot(212)
			pl.plot(time_a,airmass,'r.')
			
			# get the template location
			imx,imy,tx,ty=GetPos(starid_s[j],section)
			print "star %d [%d:%.2f,%.2f]" % ((j+1),starid_s[j],tx,ty)
			
			pl.show()
	
			print len(time_t[j])

	return exp_time,rms_med_final,bins_longer,binsr,flux_new,rms_med_final_n,binsr_n


def CustomRMSvsBin(exptime):
	
	# 2011-08-04
	#files=['/Volumes/DATA/NITES/data/2011-08-04/reduced/phot/A6/1_1/star_00169_1_1_2011-08-04.lc.txt','/Volumes/DATA/NITES/data/2011-08-04/reduced/phot/A6/1_2/star_00856_1_2_2011-08-04.lc.txt','/Volumes/DATA/NITES/data/2011-08-04/reduced/phot/A6/1_2/star_00682_1_2_2011-08-04.lc.txt','/Volumes/DATA/NITES/data/2011-08-04/reduced/phot/A6/2_2/star_00454_2_2_2011-08-04.lc.txt','/Volumes/DATA/NITES/data/2011-08-04/reduced/phot/A6/2_2/star_00640_2_2_2011-08-04.lc.txt','/Volumes/DATA/NITES/data/2011-08-04/reduced/phot/A6/2_2/star_00101_2_2_2011-08-04.lc.txt','/Volumes/DATA/NITES/data/2011-08-04/reduced/phot/A6/2_2/star_00696_2_2_2011-08-04.lc.txt','/Volumes/DATA/NITES/data/2011-08-04/reduced/phot/A6/2_2/star_00452_2_2_2011-08-04.lc.txt','/Volumes/DATA/NITES/data/2011-08-04/reduced/phot/A6/2_2/star_00945_2_2_2011-08-04.lc.txt']
	
	# 2011-08-12
	files=['/Volumes/DATA/NITES/data/2011-08-12/reduced/phot/A6/1_1/star_00169_1_1_2011-08-12.lc.txt','/Volumes/DATA/NITES/data/2011-08-12/reduced/phot/A6/1_2/star_00856_1_2_2011-08-12.lc.txt','/Volumes/DATA/NITES/data/2011-08-12/reduced/phot/A6/1_2/star_00682_1_2_2011-08-12.lc.txt','/Volumes/DATA/NITES/data/2011-08-12/reduced/phot/A6/2_2/star_00454_2_2_2011-08-12.lc.txt','/Volumes/DATA/NITES/data/2011-08-12/reduced/phot/A6/2_2/star_00640_2_2_2011-08-12.lc.txt','/Volumes/DATA/NITES/data/2011-08-12/reduced/phot/A6/2_2/star_00101_2_2_2011-08-12.lc.txt','/Volumes/DATA/NITES/data/2011-08-12/reduced/phot/A6/2_2/star_00696_2_2_2011-08-12.lc.txt','/Volumes/DATA/NITES/data/2011-08-12/reduced/phot/A6/2_2/star_00452_2_2_2011-08-12.lc.txt','/Volumes/DATA/NITES/data/2011-08-12/reduced/phot/A6/2_2/star_00945_2_2_2011-08-12.lc.txt']
	
	for i in range(0,len(files)):
		
		if i == 0:
			time_t,flux_t=ReadStar(files[i])
			lim=len(time_t)
		if i > 0:
			time,flux=ReadStar(files[i])
			time_t=np.vstack([time_t,time[:lim]])
			flux_t=np.vstack([flux_t,flux[:lim]])
	
	# define the bins for each star
	# bin until there are only 5 binned points
	z=1
	bins=[0]*(len(time_t[0])/5)
	for i in range(0,len(bins)):
		bins[i]=z
		z=z+1
	
	# un-normalised data
	rms_d_b=np.empty((len(time_t),len(bins)))
	exp_time=np.empty((len(time_t),len(bins)))
			
	for k in range(0,len(time_t)):
		
		for j in range(0,len(bins)):
			d_b=np.empty(len(time_t[k])/bins[j])
		
			q=0
			for i in range(0,len(time_t[k])/bins[j]):
				d_b[i]=((np.sum(flux_t[k][q:q+bins[j]]))/bins[j])
				q=q+bins[j]
		
			# rms in mag
			rms_d_b[k][j]=(np.std(d_b)/np.average(flux_t[k]))
			# exp time in secs per bin
			exp_time[k][j]=bins[j]*exptime
	
	# now median out the rms vs bin size from the 9 stars
	rms_med=np.copy(rms_d_b.transpose())
	rms_med_final=np.median(rms_med,axis=1)
	
	# now workout the 1/N^0.5 line for plotting gausian noise
	bins_longer=np.arange(1,201,1)
	binsr=(1/np.sqrt(bins_longer))*rms_med_final[0]
	
	pl.figure(40)
	pl.title('RMS vs Exposure Time')
	pl.ylabel('RMS')
	pl.xlabel('Exposure Time (s)')
	pl.loglog(exp_time[0],rms_med_final,'r.',(bins_longer)*exptime,binsr,'k--')
	pl.show()
	
	from datetime import date
	
	date=date.today()	
	name='RMSvsBin_TOTAL_%s.lc.txt' % (date)
	name2='RMSvsBin_1n_TOTAL_%s.lc.txt' % (date)

	z=np.concatenate((exp_time[0],rms_med_final)).reshape(2,len(exp_time[0])).transpose()
	np.savetxt(name,z,fmt='%.3f    %.5f')
	
	z2=np.concatenate(((bins_longer)*exptime,binsr)).reshape(2,len(binsr)).transpose()
	np.savetxt(name2,z2,fmt='%.3f    %.5f')
	
	
	return time_t,flux_t


def GetPos(starid,section):
	
	command='/Volumes/DATA/NITES/Results/M71/Photometry/Batch1/A8.0D12S8C12/%s/M71_%s.coo' % (section, section)
	
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
			
			print "X Y: %.2f %.2f" % (xk[0],yk[0])
	
	# get image subsection coords into tpl full 
	# frame coords for RA and DEC calculations
	if section == '1_1':
		# x(tplr1_1) = x(tpl) - 20 ; y(tplr1_1) = y(tpl) - 20
		xk[0]=xk[0]-20.0
		yk[0]=yk[0]-20.0
		print "x(tp1): %.2f y(tpl): %.2f" % (xk[0],yk[0])
	
	if section == '1_2':
		# x(tplr1_2) = x(tpl) - 20 ; y(tplr1_2) = y(tpl) + 492
		xk[0]=xk[0]-20.0
		yk[0]=yk[0]+492.0
		print "x(tp1): %.2f y(tpl): %.2f" % (xk[0],yk[0])

	if section == '2_1':
		# x(tplr1_2) = x(tpl) + 492 ; y(tplr1_2) = y(tpl) - 20
		xk[0]=xk[0]+492.0
		yk[0]=yk[0]-20.0
		print "x(tp1): %.2f y(tpl): %.2f" % (xk[0],yk[0])
	
	if section == '2_2':
		# x(tplr1_2) = x(tpl) + 492 ; y(tplr1_2) = y(tpl) + 492
		xk[0]=xk[0]+492.0
		yk[0]=yk[0]+492.0
		print "x(tp1): %.2f y(tpl): %.2f" % (xk[0],yk[0])
	
	return imx,imy,xk[0],yk[0]


def OutputRMSvsFluxFile(section,flux,rms):

	from datetime import date
	
	date=date.today()	
	name='RMSvsFlux_%s_%s.lc.txt' % (section,date)

	z=np.concatenate((flux,rms)).reshape(2,len(flux)).transpose()
	np.savetxt(name,z,fmt='%.3f    %.5f')

	return 0

def OutputRMSvsMagFile(section,mag,rms):

	from datetime import date
	
	date=date.today()	
	name='RMSvsMag_%s_%s.lc.txt' % (section,date)

	z=np.concatenate((mag,rms)).reshape(2,len(mag)).transpose()
	np.savetxt(name,z,fmt='%.3f    %.5f')

	return 0

def OutputNoiseModelFile(fluxes, object_model, sbg_model, rn_model, dc_model, ff_model, scint_model, total_err_model_wscintff):
	
	from datetime import date
	
	date=date.today()	
	name='NoiseModel_Total_%s.lc.txt' % (date)
	
	z=np.concatenate((fluxes,object_model,sbg_model,rn_model,dc_model,ff_model,scint_model,total_err_model_wscintff)).reshape(8,len(fluxes)).transpose()
	np.savetxt(name,z,fmt='%.3f    %.5f    %.5f    %.5f    %.5f    %.5f    %.5f    %.5f')

	return 0



########
# Main #
########

exptime=30.0
gain=1.22
section='1_1'

# get a list of the images
templist=cmd.getoutput('ls /Volumes/DATA/NITES/data/2011-08-04/reduced/M71-*.fits').split('\n')

# get sky background counts
trim=0
sky_bg_pp=GetSkyMedian(gain,trim) 
#sky_bg_pp=100.0 # to make model fit better 20120804 = 100
#sky_bg_pp=1 # Dark 20110804 181.4 Bright 20110812 805.3 

# get airmasses
airmass=np.empty(len(templist))
time_a=np.empty(len(templist))

for i in range(0,len(templist)):
	hdulist = pf.open(templist[i])
	airmass[i]=hdulist[0].header['AIRMASS']
	time_a[i]=hdulist[0].header['HJD-MID']
	
fluxes, sbg_model, object_model, rn_model, dc_model, ff_model, scint_model, total_err_model_wff, total_err_model_wscintff = GetNoiseModel(airmass,exptime,sky_bg_pp,gain) 

# get template fluxes and errs
tflux,tfluxerr=GetFlux(section)

# get rms's
n,n2,rms,plot_flux,starid=GetRMSvsFlux(section,tflux)

# output rms vs flux file
#out=OutputRMSvsFluxFile(section,plot_flux[n],rms[n])

# output the noise model file
out2=OutputNoiseModelFile(fluxes, object_model, sbg_model, rn_model, dc_model, ff_model, scint_model, total_err_model_wscintff)

# get rms vs bin size
#exp_time,rms_med_final,bins_longer,binsr,flux_new,rms_med_final_n,binsr_n=GetRMSvsBin(plot_flux,rms,starid,n,exptime,time_a,airmass)

#time_tc,flux_tc=CustomRMSvsBin(exptime)

# find mag at 1% level for planet hunt
mags=(-2.5*np.log10(plot_flux[n]))+26.974

# output rms vs mag file
#out=OutputRMSvsMagFile(section,mags,rms[n])

line1=[13.0,14.0,15.0,16.0,17.0,18.0,19.0]
line2=[0.01,0.01,0.01,0.01,0.01,0.01,0.01]
line3=[16.5,16.5,16.5,16.5,16.5,16.5,16.5]
line4=[0.0001,0.001,0.005,0.01,0.05,0.1,0.5]
line5=[17.5,17.5,17.5,17.5,17.5,17.5,17.5]

# plots
pl.figure(1)
pl.loglog(fluxes, object_model, 'g-',fluxes, sbg_model, 'b-',fluxes, rn_model, 'r-', fluxes, dc_model, 'm-', fluxes, ff_model, 'c--', fluxes, scint_model, 'y-', fluxes, total_err_model_wff, 'k-', fluxes, total_err_model_wscintff, 'k--')
pl.legend(('Object', 'Sky', 'Read Noise', 'Dark Current', 'Flat Fields', 'Scintillation', 'Total_CCD+ff', 'Total_CCD+ff+Scint' ), loc=0)
pl.loglog(plot_flux[n],rms[n],'r.',plot_flux[n2],rms[n2],'r+')
pl.ylim(0.0001,1)
pl.xlim(1000,1000000)
pl.ylabel('RMS (mags)')
pl.xlabel('Flux')

pl.figure(2)
pl.title('RMS vs Mag (binned 5min)')
pl.ylabel('RMS (mag)')
pl.xlabel('V Mag')
pl.semilogy(mags,rms[n],'r.',line1,line2,'k-',line3,line4,'k-', line5,line4,'k-')
			
pl.show()
