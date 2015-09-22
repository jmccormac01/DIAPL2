# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_CycleNights.py - 	a python program to cycle through the nights. 
#
#		This is primarily to look for SX PHe type objects. Can be used to exclude 
# 		bad nights also				
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 05/12/11 - 	code writen
# 05/12/11 -	code tested
#				
#

import commands as cmd
import numpy as np
import pylab as pl
import os, os.path


def ReadStar(file):
	
	time,flux,err,sky=np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
		
	return time,flux,err,sky


def GetStartsEnds(time):
	
	starts=[0]
	ends=[]
	
	for j in range(0,len(time)):
		if j < len(time)-1:
			if time[j+1] - time[j] > 0.5:
				starts.append(j+1)
				ends.append(j)
	
	ends.append(len(time))
	
	return starts,ends


def Fit(time,flux,err):
	
	pvals=np.empty(len(time))
	flux_new=np.empty(len(time))
	mag=np.empty(len(time))

	# best line of fit for flux out of transit
	coeffs=np.polyfit(time,flux,2)
	
	# best vals for mags
	besty=np.polyval(coeffs,time)
	
	# finding p to fit the data -- equation -- p2xk^2 + p1xk + p0 =yk
	for i in range(0,len(flux)):
		pvals[i]=coeffs[2]+(coeffs[1]*time[i])+(coeffs[0]*(time[i]*time[i]))

	for i in range(0,len(flux)):
		flux_new[i]=flux[i]/pvals[i]
	
	mag=-2.5*np.log10(flux_new)
	
	mag_err=err/np.average(flux)
	
	rms=np.std(flux_new)
	
	pl.figure(20)
	pl.subplots_adjust(hspace=0.5)
	pl.subplot(211)	
	pl.title('Section Fit')
	pl.errorbar(time,flux,yerr=err,fmt='r.')
	pl.plot(time,besty,'k-')
	
	pl.subplot(212)	
	pl.title('Section Fitted')
	pl.errorbar(time,mag,yerr=mag_err,fmt='r.')
	pl.gca().invert_yaxis()
	pl.show()
	
	# CLIP!!!
	
	return rms,flux_new,mag,mag_err


def PrintFile(name,time,flux,err):
	
	savedir="variables_final"
	
	if os.path.exists(savedir) == False:
		os.mkdir(savedir)
	
	z=np.concatenate((time,flux,err)).reshape(3,len(time)).transpose()
	np.savetxt(name,z,fmt='%.8f    %.3f    %.3f')

	return 0
	
	
def IndexGnuplot(file,time,flux,err,sky,starts,ends,toggle):
	
	# create 2 blank lines between nights for gnuplot point exclusion
	if toggle == 0:
		print "Replacing %s" % (file)
		comm="rm -rf %s" % (file) 
		os.system(comm)
		

		f=open(file,'w')
		for j in range(0,len(starts)):
			i_time=time[starts[j]:ends[j]]
			i_flux=flux[starts[j]:ends[j]]
			i_err=err[starts[j]:ends[j]]
			i_sky=sky[starts[j]:ends[j]]
	
			for k in range(0,len(i_time)):
				line = "%.8f    %.3f    %.3f   %.3f\n" % (i_time[k],i_flux[k],i_err[k],i_sky[k])
				f.write(line)
			f.write('\n\n')	
	
	# write something after to clean up the blank spaces for the final lc
	if toggle == 1:
		print "Replacing %s" % (file)
		comm="rm -rf %s" % (file) 
		os.system(comm)
		
		z=np.concatenate((time,flux,err,sky)).reshape(4,len(time)).transpose()
		np.savetxt(file,z,fmt='%.8f    %.3f    %.3f    %.3f')
	
	return 0


templist=cmd.getoutput('ls star_01073*').split('\n')
section="%s_%s" % (templist[0].split('_')[2],templist[0].split('_')[3])
starid=np.empty(len(templist))

# toggles
cycle = 0
analyse = 0
index = 0
unindex = 1

for i in range(0,len(templist)):
	
	time,flux,err,sky=ReadStar(templist[i])

	starid[i]=int(templist[i].split('_')[1])
	
	# get start and end times of each night
	starts,ends=GetStartsEnds(time)
	
	if index != unindex:
		if index == 1:
			done=IndexGnuplot(templist[i],time,flux,err,sky,starts,ends,0)
			
		if unindex == 1:
			done=IndexGnuplot(templist[i],time,flux,err,sky,starts,ends,1)
		
	# cycle through the nights?
	if cycle == 1:
	
		# kept data arrays
		time_keep=np.empty(0)
		flux_keep=np.empty(0)
		err_keep=np.empty(0)
		sky_keep=np.empty(0)
		
		for j in range(0,len(starts)):
		
			jd_n=int(time[starts[j]])
			p_time=time[starts[j]:ends[j]]-jd_n
			p_flux=flux[starts[j]:ends[j]]
			p_err=err[starts[j]:ends[j]]
			p_sky=sky[starts[j]:ends[j]]
			
			# bin data into 1 min bins
			#bin=len(p_flux)/3
			#bf=int(len(p_time)/bin)
			# binned arrays
			#tbinned=np.empty(bin)
			#fbinned=np.empty(bin)
			#errbinned=np.empty(bin)
			#	
			#for k in range(0,len(tbinned)):
			#	tbinned[k]=np.average(p_time[(k*bf):bf*(k+1)])
			#	fbinned[k]=np.average(p_flux[(k*bf):bf*(k+1)])
			#	errbinned[k]=np.average(p_err[(k*bf):bf*(k+1)])
			
			print "%d\t\t%d" % ((j+1),jd_n)
			
			pl.figure(20)
			pl.errorbar(p_time,p_flux,yerr=p_err,fmt='r.')
			pl.title("Star: %d" % (starid[i]))
			pl.ylabel('Differential Flux')
			pl.xlabel('JD %d+' % (jd_n))	
			pl.show()
			
			done=0
			if analyse > 0:
				while done != 1:
					q=raw_input("Action: ")
					
					# fit a night and get rms
					if str(q) == 'f':
						rms,flux_new,mag,mag_err=Fit(p_time,p_flux,p_err)
						print "RMS %.6f" % (rms)
						done = 1
						
					# keep a full night	
					if str(q) == 'k':
						# keep the real times and not the trimmed ones
						time_keep=np.hstack((time_keep,time[starts[j]:ends[j]]))
						flux_keep=np.hstack((flux_keep,p_flux))
						err_keep=np.hstack((err_keep,p_err))
						sky_keep=np.hstack((sky_keep,p_sky))
						done = 1
						
					# do nothing
					if str(q) == "":
						done = 1
						continue
		
			if j == len(starts)-1:
				print "Output final file?"
				
				# flux - normal data
				yn=raw_input("(e.g. y): ")
				if str(yn) == 'y':
					name="star_%05d_%s_wEcS_f2.lc.txt" % (starid,section)
					out=PrintFile(name,time_keep,flux_keep,err_keep)

				
				
		
	