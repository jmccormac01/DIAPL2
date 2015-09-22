
# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_NormalizeSXPhe.py - 	a program to normalize SXPhe variations
#										
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 09/12/11 - 	code writen
# 09/12/11 -	code tested
#

import commands as cmd
import numpy as np
import pylab as pl

def ReadStar(file):
	
	time,flux,err=np.loadtxt(file,usecols=[0,1,2],unpack=True)
		
	return time,flux,err

def Fit(time,flux,err):
	
	pvals=np.empty(len(time))
	flux_new=np.empty(len(time))
	mag=np.empty(len(time))

	# best line of fit for flux out of transit
	coeffs=np.polyfit(time,flux,1)
	
	# best vals for mags
	besty=np.polyval(coeffs,time)
	
	# finding p to fit the data -- equation -- p2xk^2 + p1xk + p0 =yk
	for i in range(0,len(flux)):
		pvals[i]=coeffs[1]+(coeffs[0]*time[i])

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
	
	return rms,flux_new,mag,mag_err

	
def PrintFile(name,time,flux,err,mag,mag_err,toggle):
	
	if toggle == 0:
		z=np.concatenate((time,flux,err,mag,mag_err)).reshape(5,len(time)).transpose()
		np.savetxt(name,z,fmt='%.8f    %.3f    %.3f    %.3f    %.3f')
	if toggle == 1:
		z=np.concatenate((time,mag,mag_err)).reshape(3,len(time)).transpose()
		np.savetxt(name,z,fmt='%.8f    %.3f    %.3f')

	return 0
	
	
templist=cmd.getoutput('ls SXPhe_*').split('\n')

total_time=np.empty(0)
total_flux=np.empty(0)
total_err=np.empty(0)
total_mag=np.empty(0)
total_mag_err=np.empty(0)


for i in range(0,len(templist)):
	time,flux,err=ReadStar(templist[i])
	rms,flux_new,mag,mag_err=Fit(time,flux,err)
	
	total_time=np.hstack((total_time,time))
	total_flux=np.hstack((total_flux,flux))
	total_err=np.hstack((total_err,err))
	total_mag=np.hstack((total_mag,mag))
	total_mag_err=np.hstack((total_mag_err,mag_err))
	
	out=PrintFile(templist[i],time,flux,err,mag,mag_err,0)

out=PrintFile('SXPhe_108_total.lc.txt',total_time,total_flux,total_err,total_mag,total_mag_err,1)
	
		