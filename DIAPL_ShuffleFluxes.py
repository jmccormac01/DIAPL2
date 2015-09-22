# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_ShuffleFluxes.py - 	a program to shuffle the flux values for FAP
#										
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 21/11/11 - 	code writen
# 21/11/11 -	code tested

import commands as cmd
import random
import numpy as np


def GetFlux(section):
	
	command='/Volumes/DATA/NITES/Results/M71/Photometry/Batch1/A8.0D12S8C12/%s/M71_%s.flux' % (section, section)
	
	f=open(command).readlines()
	
	tflux=np.empty(len(f))
	tfluxerr=np.empty(len(f))
	
	for i in range(0,len(f)):
		tflux[i]=float(f[i].split('\n')[0].split()[3])
		tfluxerr[i]=float(f[i].split('\n')[0].split()[4])


	return tflux,tfluxerr
	
	
templist=cmd.getoutput('ls star_*.lc.txt').split('\n')
section="%s_%s" % (templist[0].split('_')[2],templist[0].split('_')[3][0])

# get the sextractor fluxes
tflux,tfluxerr=GetFlux(section)


for i in range(0,len(templist)):
	
	starid=int(templist[i].split('_')[1])
	
	time=np.loadtxt(templist[i], usecols=[0], unpack=True)
	
	name="%s_shuf.lc.txt" % (templist[i].split('.')[0])
	
	flux_new=np.empty(len(time))
	for j in range(0,len(time)):
		flux_new[j]=tflux[(starid-1)] + (tfluxerr[(starid-1)]*random.gauss(0,1))
	
	
	z=np.concatenate((time,flux_new)).reshape(2,len(time)).transpose()
	np.savetxt(name,z,fmt='%.8f    %.3f')
	
	print "[Shuffle] %d/%d" % ((i+1),len(templist))
	
	