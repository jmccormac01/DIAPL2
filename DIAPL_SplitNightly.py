

# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_SplitNightly.py - 	a program to split the long star_* files into a single 
#							night
# 						
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 02/12/11 - 	code writen
# 02/12/11 -	code tested
#

import commands as cmd
import numpy as np
import pylab as pl

def ReadStar(file):

	time,flux=np.loadtxt(file,usecols=[0,1],unpack=True)
		
	return time,flux


def GetFlux(section):
	
	file='/Volumes/DATA/NITES/Results/M71/Photometry/Batch1/A8.0D12S8C12/%s/M71_%s.flux' % (section, section)
	tflux,tfluxerr=np.loadtxt(file,usecols=[3,4],unpack=True)

	return tflux,tfluxerr


def OutputFile(starid,section,time,flux,note):
	
	name='star_%05d_%s_f1_b_%s.lc.txt' % (starid,section,note)

	z=np.concatenate((time,flux)).reshape(2,len(flux)).transpose()
	np.savetxt(name,z,fmt='%.8f    %.3f')

	return 0


########
# MAIN #
########

section = '1_1'
note='2011-08-04'

run_split = 1
check_split = 1

if run_split > 0:
	comm='ls /Volumes/DATA/NITES/Results/M71/Photometry/Batch1/BLS/%s/binned/star_*.lc.txt' % (section)
	
	templist=cmd.getoutput(comm).split('\n')
	
	# get template fluxes and errs
	tflux,tfluxerr=GetFlux(section)
	
	starid=np.empty(len(templist))
		
	# get the stddev for all
	for i in range(0,len(templist)):
		time,flux=ReadStar(templist[i])
	
		# add the template flux from Sextractor	
		starid[i]=int(templist[i].split('_')[2])
		flux=flux+tflux[(starid[i]-1)]
			
		# get data on 2011-08-04
		# 2455778.39028
		# 2455778.67986		
		# get data on 2011-08-12
		# 2455786.45356
		# 2455786.64055
		
		c1=np.empty(len(time))
		c2=np.empty(len(time))
		for j in range(0,len(time)):
			c1[j]=abs(time[j]-2455778.39028)
			c2[j]=abs(time[j]-2455778.67986)
	
		start=np.where(c1==min(c1))[0][0]
		end=np.where(c2==min(c2))[0][0]
		
		# output the data from chosen night
		out=OutputFile(starid[i],section,time[start:end],flux[start:end],note)
		
		print "[2011-08-04] %d/%d" % ((i+1),len(templist))

# code something to check the new files lengths so they are all roughly the same
if check_split > 0:
	comm2='ls star_*_%s.lc.txt' % (note)
	t=cmd.getoutput(comm2).split('\n')

	f=np.empty(len(t))
	for i in range(0,len(t)):
		f[i]=len(open(t[i]).readlines())
	
	pl.plot(f,'r-')
	pl.show()


