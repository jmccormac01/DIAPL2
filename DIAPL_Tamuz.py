
# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_Tamuz.py - 	a program to detrend DIAPL2 reduced data
#										
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 16/11/11 - 	code writen
# 16/11/11 -	code tested
#
#		extract lc's 
# 		get the corresponding airmasses - uses 1's for now
#			



from numpy import * 
import numpy as np
from numarray import *                       
import pyfits, os, sys, string
from math import *
from math import sqrt
import time
import commands
from pylab import *
import pylab
import matplotlib.pyplot as plt
from matplotlib import axes 
from datetime import date
import random


def GetFlux(section):
	
	command='/Volumes/DATA/NITES/Results/M71/Photometry/Batch1/A8.0D12S8C12/%s/M71_%s.flux' % (section, section)
	
	f=open(command).readlines()
	
	tflux=np.empty(len(f))
	tfluxerr=np.empty(len(f))
	
	for i in range(0,len(f)):
		tflux[i]=float(f[i].split('\n')[0].split()[3])
		tfluxerr[i]=float(f[i].split('\n')[0].split()[4])


	return tflux,tfluxerr

def beta_min(ff,AIRMASS,nim):

	rij_new = None
	ffn = None

	# original data structures
	rms_ffo=np.empty(len(ff))
	av_ffo=np.empty(len(ff))

	rms_ff=np.empty(len(ff))
	rms_ffn=np.empty(len(ff))

	# shuffle
	#ff_shuffle=np.empty( (len(ff), len(templist) ) )

	# arrays for c
	c=zeros(len(ff),float)
	sum_top_c=zeros(nim,float)
	sum_bottom_c=zeros(nim,float)
	
	# arrays for a
	a=zeros(nim,float)
	sum_top_a=zeros(len(ff),float)
	sum_bottom_a=zeros(len(ff),float)
	
	# arrays for s
	Ssquare=zeros(len(ff), float)
	Ssquare_j=zeros(nim,float)

	t=zeros(nim,float)
	s=zeros(len(ff),float)

	#chi sq arrays
	chisq=zeros(1,float)
	
	ff_shuffle=np.copy(ff)
	#for i in range(0,len(ff)):
	#	ff_shuffle[i]=ff[i]
	
	for i in range(0,len(ff_shuffle)):
		random.shuffle(ff_shuffle[i])
		
	print "Average ff_shuffle[0]: " + str(average(ff_shuffle[0]))
	
	print "\nBeta_min: Gathering uncertanties for " + str(len(ff)*nim) + " points"
	
	sigma = np.sqrt(np.abs(ff_shuffle))
	sigmasq = sigma ** 2
	
	#sigma=np.empty( (len(ff), nim) )
	#for i in range(0, len(ff)):
	#	for j in range(0, len(templist)):
	#		sigma[i,j]=sqrt(ff_shuffle[i,j])		

	for i in range(0, len(ff_shuffle)):
		av_ffo[i]=average(ff_shuffle[i])
		rms_ffo[i]=std(ff_shuffle[i])/average(ff_shuffle[i])	


	for trend in range(0, 1):

		OLD_SS = 1e38
		DIFF = 1
		
		counter = 1
		
		# test this when the code has been updated
		a=AIRMASS
		#a=ones(len(templist),float)
		
		print "\nDetrending trend number " + str(trend+1) + '\n'
		
		while DIFF > 0.0001:
			print "Beginning Iteration Number: " + str(counter)
			
			print "Calculating c"
			c = ((rij * a)/sigmasq).sum(axis = 1) / ((a**2) / sigmasq).sum(axis = 1)	
			#for i in range (0, len(ff)):
			#	for j in range(0, len(templist)):
			#		sum_top_c[j]=rij[i,j]*a[j]/(pow(sigma[i,j],2))
			#		sum_bottom_c[j]=(pow(a[j],2))/(pow(sigma[i,j],2))
			#	c[i]=sum(sum_top_c)/sum(sum_bottom_c) # Equation (2) from Tamuz Paper
			
			print "Calculating new a"
			ccol = c[:,np.newaxis]
			a = ((rij * ccol) / sigmasq).sum(axis = 0) / (ccol ** 2 / sigmasq).sum(axis = 0)
			# calculate the new effective airmass for each frame by
			#for j in range(0, len(templist)):
			#	for i in range(0, len(ff)):
			#		sum_top_a[i]=rij[i,j]*c[i]/(pow(sigma[i,j],2))
			#		sum_bottom_a[i]=(pow(c[i],2))/(pow(sigma[i,j],2))
			#	a[j]=sum(sum_top_a)/sum(sum_bottom_a)	# Equation (4) from Tamuz Paper a(1)_j	
		
			
			cta = c[:,np.newaxis] * a
			Ssquare=(((rij-cta)**2)/sigmasq).sum(axis = 1)
			# now make an estimate of S^2
			#for i in range(0, len(ff)):
			#	for j in range(0, len(templist)):
			#		Ssquare_j[j]=pow((rij[i,j] - c[i]*a[j]),2)/(pow((sigma[i,j]),2)) 
			#
			#	Ssquare[i]=sum(Ssquare_j)
		
			Ss=sum(Ssquare)
		
			print str(OLD_SS) + ' ' + str(Ss)
		    
			DIFF = fabs((OLD_SS - Ss)/OLD_SS)
			
			print "Difference is : " + str(DIFF)     
			counter = counter + 1
			print "Next Iteration Number : " + str(counter)
		
			OLD_SS = Ss
			
			print 'Ss = ' + str(Ss) + '\n'
			
			if DIFF<0.0001:
				chisq[trend]=OLD_SS
				print "Trend " + str(trend) + " final chi squared: " + str(OLD_SS) + " [" + str(chisq[trend]) + "] \n"
			
		#new residuals by subtracting the linear effects found above
		for i in range(0, len(ff)):
			if rij_new is None:
			    rij_new = np.empty( (len(ff), nim ) )
			for j in range(0, nim):
				rij_new[i,j] = rij[i,j] - (c[i]*a[j])
		
		
		#list of new flux values 		
		for i in range(0, len(ff)):
			if ffn is None:
			    ffn = np.empty( (len(ff), nim ) )
			
			ffn[i]=(rij_new[i])
			
	
		for i in range(0, len(ff)):
			rms_ffn[i]=std(ffn[i])/average(ff[i])
	
		
		beta_shuffle=(rms_ffo-rms_ffn)/rms_ffo

		beta_sort=sort(beta_shuffle)
	
		# fraction where we want to chose beta_min
		alpha = 0.90
		
		beta_min_loc=int(len(beta_sort)*alpha)
		
		figure(0)
		title('Probability vs Beta')
		pylab.plot(beta_sort,'r.')
		xlabel('Number')
		ylabel('Beta')
		
		betamin=beta_sort[beta_min_loc]
		
		print 'Beta_min: ' + str(betamin) 
		print 'Beta_min @ ' + str(beta_min_loc)
	
	return betamin



# len(templist) == no. stars == len(ff)
# nim == no. images

# get image list
templist=commands.getoutput('ls star_*.lc.txt').split('\n')
# get image subsection
section=str(templist[0].split('.')[0][-3:])
# get the sextractor fluxes
tflux,tfluxerr=GetFlux(section)

gain=1.22
nim=25741 # 677 images on 2011-08-04 1_1

starid=np.empty(len(templist))
flux=np.empty( (len(templist), nim) )
hjd=np.empty( (len(templist), nim) )

# get fluxes ready
for i in range(0,len(templist)):
	starid[i]=int(templist[i].split('_')[1])
	f=open(templist[i]).readlines()
	flux[i] = np.loadtxt(templist[i], usecols=[1]) + tflux[starid[i] - 1]
	hjd[i] = np.loadtxt(templist[i], usecols=[0])

	# [15752:16429] - was figure out from a simple analysis to show only data from 2011-08-04 for 1_1
	# I checked that all lc's have the same HJD vales for this range manually 
	

	print "[Get Stars] %d/%d..." % ((i+1),len(templist))

# scale by the gain
flux=flux*gain

bad=np.array([0,3,18,32,48,49,70,82,118,121,122,126,127,150,177,198,227,233,259,293,304,330,336,344,345,355,377,378,397,399,435,452,456,466,469,471,474,476,492,496,505,519,550,572,585,629,634,676,701,714,716,854]) # 2011-08-04

zeros_1=np.zeros(len(flux))
zeros_1[bad]=1.0
index=np.where(zeros_1==0.)

ff=flux[index]

# get errors
print "Gathering uncertanties for " + str(len(ff)*nim) + " points"
#sigma=np.empty( (len(ff), nim) )
#for i in range(0, len(ff)):
#	for j in range(0, nim):
#		sigma[i,j]=sqrt(abs(ff[i,j]))	
#	
#	print "[Get Errors] %d/%d..." % ((i+1),len(ff))

sigma = np.sqrt(np.abs(ff))

rij_new = None
ffn = None

# original data structures
rms_ffo=np.empty(len(ff))
av_ffo=np.empty(len(ff))

rms_ff=np.empty(len(ff))
rms_ffn=np.empty(len(ff))

# arrays for c
c=zeros(len(ff),float)
sum_top_c=zeros(nim,float)
sum_bottom_c=zeros(nim,float)
	
# arrays for a
a=zeros(nim,float)
sum_top_a=zeros(len(ff),float)
sum_bottom_a=zeros(len(ff),float)
	
# arrays for s
Ssquare=zeros(len(ff), float)
Ssquare_j=zeros(nim,float)

t=zeros(nim,float)
s=zeros(len(ff),float)

#chi sq arrays
chisq=zeros(10,float)
	
stddev=np.empty(len(ff))
stddev_o=np.empty(len(ff))
rms=np.empty(len(ff))
rms_o=np.empty(len(ff))
plot_flux=np.empty(len(ff))


#################################################################################################
					# Begin the looping over iterations from here #
#################################################################################################

sigmasq = sigma ** 2
rij = np.copy(ff)
#rij = None
#for i in range(0, len(ff)):
	#if rij is None:
	#	    rij = np.empty( ff.shape )

	# rij done for each star, rows are per star, columns are per image
	# multiply by the gain so we are working in electrons
	#rij[i]=ff[i]

a=ones(nim,float)

#b_min=beta_min(ff,a,nim)


run_t = 1

if run_t > 0:

	trend = 0
	stop = 10
	#while stop > 1:
	
	for trend in range(0,1):
		OLD_SS = 1e38
		DIFF = 1
		
		counter = 1
		trend = trend + 1
		
		# test this when the code has been updated
		#a=AIRMASS
		
		print "\nDetrending trend number " + str(trend) + '\n'
		
		while DIFF > 0.0001:
			print "Beginning Iteration Number: " + str(counter)
			
			print "Calculating c"
			c = ((rij * a)/sigmasq).sum(axis = 1) / ((a**2) / sigmasq).sum(axis = 1)
			#for i in range (0, len(ff)):
			#	for j in range(0, nim):
			#		sum_top_c[j]=rij[i,j]*a[j]/(sigmasq[i,j])
			#		sum_bottom_c[j]=(a[j]*a[j])/(sigmasq[i,j])
			#	c[i]=sum(sum_top_c)/sum(sum_bottom_c) # Equation (2) from Tamuz Paper
			#	
			#	print "[c] %d/%d: %f" % ((i+1),len(ff),c[i])  	
			
			print "Calculating new a"
			ccol = c[:,np.newaxis]
			a = ((rij * ccol) / sigmasq).sum(axis = 0) / (ccol ** 2 / sigmasq).sum(axis = 0)
			# calculate the new effective airmass for each frame by
			#for j in range(0, nim):
			#	for i in range(0, len(ff)):
			#		sum_top_a[i]=rij[i,j]*c[i]/(sigmasq[i,j])
			#		sum_bottom_a[i]=(c[i]*c[i])/(sigmasq[i,j])
			#	a[j]=sum(sum_top_a)/sum(sum_bottom_a)	# Equation (4) from Tamuz Paper a(1)_j	
			#	
			#	print "[a] %d/%d: %f" % ((j+1),nim,a[j])
		
			# now make an estimate of S^2
			
			cta = c[:,np.newaxis] * a
			Ssquare=(((rij-cta)**2)/sigmasq).sum(axis = 1)
			
			#for i in range(0, len(ff)):
			#	for j in range(0, nim):
			#		Ssquare_j[j]=((rij[i,j] - (c[i]*a[j]))*(rij[i,j] - (c[i]*a[j])))/(sigmasq[i,j]) 
			#
			#	Ssquare[i]=sum(Ssquare_j)
			#	
			#	print "[S^2] %d/%d" % ((i+1),len(ff))
				
			Ss=sum(Ssquare)
		
			print str(OLD_SS) + ' ' + str(Ss)
		    
			DIFF = fabs((OLD_SS - Ss)/OLD_SS)
			
			print "Difference is : " + str(DIFF)     
			counter = counter + 1
			print "Next Iteration Number : " + str(counter)
		
			OLD_SS = Ss
			
			print 'Ss = ' + str(Ss) + '\n'
			
			if DIFF<0.0001:
				chisq[trend]=OLD_SS
				print "Trend " + str(trend) + " final chi squared: " + str(OLD_SS) + " [" + str(chisq[trend]) + "] \n"
		
		
		cta_f = c[:,np.newaxis] * a	
		rij_new=rij-cta_f
		
		#new residuals by subtracting the linear effects found above
		#for i in range(0, len(ff)):
		#	if rij_new is None:
		#	    rij_new = np.empty( (len(ff), nim ) )
		#	for j in range(0, nim):
		#		rij_new[i,j] = rij[i,j] - (c[i]*a[j])
		
		ffn=np.copy(rij_new)
		
		#list of new flux values 		
		#for i in range(0, len(ff)):
		#	if ffn is None:
		#	    ffn = np.empty( (len(ff), nim ) )
		#	
		#	ffn[i]=(rij_new[i])
			
	
		for i in range(0, len(ff)):
			rms_ffn[i]=std(ffn[i])/average(ff[i])
	
	
		if trend == 1:
			beta=(rms_ffo-rms_ffn)/rms_ffo
		if trend > 1:
			beta=(rms_ff-rms_ffn)/rms_ff
		
		for i in range(0,len(ff)):
			rms_ff[i]=rms_ffn[i]
		
		# set rij to the new residuals with the previous trend removed
		rij=rij_new
		
		clock = 0.0
		for i in range(0,len(beta)):
			if beta[i]-b_min < 0:
				clock=clock + 1.0
		
		stop_ratio = clock / float(len(beta))
		print "Stopping ratio: " + str(stop_ratio)
		
		if stop_ratio > 0.95:
			stop = 0
		
			
		# create a multiplot comparing each iteration with the original data	
		# turn on interactive plotting
		#pylab.ion()
		#figure(1)
		#if trend == 1:
		#	suptitle('RMS vs MAG')
		#	ylabel('RMS')
		#	xlabel('Flux')
		#	pylab.semilogx(av_flux_ffo,rms_ffo,'y.')
		#	pylab.semilogx(av_flux_ffo,rms_ff,'r.')
		#	ylim(0,0.1)
		#
		#if trend == 2:
		#	pylab.semilogx(av_flux_ffo,rms_ff,'k.')
		#	ylim(0,0.1)
		#	
		#if trend == 3:
		#	pylab.semilogx(av_flux_ffo,rms_ff,'b.')
		#	ylim(0,0.1)
		#	
		#if trend == 4:
		#	pylab.semilogx(av_flux_ffo,rms_ff,'g.')
		#	ylim(0,0.1)
		#	
		#figure(2)
		#if trend == 1:
		#	suptitle('Beta vs MAG')
		#	ylabel('Beta')
		#	xlabel('Flux')
		#	pylab.semilogx(av_flux_ffo,beta,'r.')
		#	ylim(0,0.1)
		#
		#if trend == 2:
		#	pylab.semilogx(av_flux_ffo,beta,'k.')
		#	
		#if trend == 3:
		#	pylab.semilogx(av_flux_ffo,beta,'b.')
		#	
		#if trend == 4:
		#	pylab.semilogx(av_flux_ffo,beta,'g.')
	
	#pylab.ioff()


for j in range(0,len(ff)):
	stddev[j]=std(ffn[j])

for j in range(0,len(ff)):
	stddev_o[j]=std(ff[j])

for i in range(0,len(ff)):
	plot_flux[i]=tflux[(starid[i]-1)]

for j in range(0,len(ff)):
	rms[j]=stddev[j]/(plot_flux[j])

for j in range(0,len(ff)):
	rms_o[j]=stddev_o[j]/(plot_flux[j])

figure(1)
semilogx(plot_flux,rms_o,'r.',plot_flux,rms,'g.')

show()

#################################################################################################
									# EYEBALL EACH LIGHTCURVE #
#################################################################################################

print "Go through each plot?"
plots_yn = raw_input("(e.g. y): ")
if str(plots_yn) == 'y':
	
	for i in range(0,len(ff)):
		print "LC: " + str(i) 
		plot(ffn[i],'g.',ff[i],'r.')
		xlim(-20,(nim+20))
		show()

print "Output files?"
files_yn = raw_input("(e.g. y): ")
if str(files_yn) == 'y':

	for i in range(0,len(templist)):
		name="star_%05d_1_1_tam.lc.txt" % (starid[i])
		
		z=np.concatenate((hjd[i],(ffn[i]+tflux[(starid[i]-1)]))).reshape(2,len(ffn[i])).transpose()
		np.savetxt(name,z,fmt='%.8f    %.3f')











