# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_CorrectErrors.py - 	a program to correct DIAPL2 errors using the 
#							SExtractor values			
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 08/12/11 - 	code writen
# 08/12/11 -	code tested
#
#

import numpy as np
import commands as cmd

###############################################
################# FUNCTIONS ###################
###############################################

def ReadStar(file):
	
	time,flux,err,sky=np.loadtxt(file,usecols=[0,1,2,3],unpack=True)
		
	return time,flux,err,sky
	
	
def GetFlux(section):
	
	command='/Volumes/DATA/NITES/Results/M71/Photometry/Batch1/A8.0D12S8C12/%s/M71_%s.flux' % (section, section)
	
	tflux,tfluxerr=np.loadtxt(command,usecols=[3,4],unpack=True)

	return tflux,tfluxerr


def GetRealErrors(err,t_err):
	
	err_new=(err/np.average(err))*t_err
	
	return err_new	
	
	
def PrintFile(name,time,flux,err,sky):
	
	z=np.concatenate((time,flux,err,sky)).reshape(4,len(time)).transpose()
	np.savetxt(name,z,fmt='%.8f    %.3f    %.3f    %.3f')

	return 0
	
	
###############################################
################### MAIN ######################
###############################################
	
# get the list of files	
templist=cmd.getoutput('ls star_*').split('\n')

# get image subsection
section="%s_%s" % (templist[0].split('_')[2],templist[0].split('_')[3][0])

# get template fluxes and errors
tflux,tfluxerr=GetFlux(section)
	
# correct all the files	
for i in range(0,len(templist)):	
	time,flux,err,sky=ReadStar(templist[i])	
	
	starid=int(templist[i].split('_')[1])
	
	err_new=GetRealErrors(err,tfluxerr[(starid-1)])
	
	name="star_%05d_%s_wEcS_f1.lc.txt" % (starid,section)
	
	out=PrintFile(name,time,flux,err_new,sky)
	print "[Corrected] %d/%d" % ((i+1),len(templist))
	
	
	
	
	
	
	