
# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_FilterStarList.py - a program to filter out bad frames using file made
#							by DIAPL_ListBadFrames.py from starbase
#
#							'batchX-badtimes.lc.txt'
#
#	RUN THIS ON R23 TO MAKE IT FASTER!!!
#
# Output:
#	star_id_section_f1.lc.txt - f1 shows filtered 1. 
#	
#	In either case of:
#		obs > 40 deg elevation and FWHM < 6.5 pixels
#		obs > 35 deg elevation and FWHM < 6.5 pixels						
#
# 	This is to stop any possible vignetting of the telescope by the dome
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 14/11/11 - 	code writen
# 14/11/11 -	code tested
# 06/12/11 - 	updated functions with np.loadtxt()
#				made program more modular
#	
#

import numpy as np
import commands as cmd
import os, os.path, sys

############################################################
###################### Functions ###########################
############################################################

def GetBadTimes():
	
	badfile='/data/waspr1/jmcc/M71-2/batch1/raw/batch1-badtimes_6.5_40.lc.txt'
	
	f=open(badfile,'r').readlines()
	
	btimes=np.empty(len(f))
	for i in range(0,len(f)):
		btimes[i]=float(f[i].split('\n')[0])
	
	#btimes=np.loadtxt(badfile,usecols=[0])
		
	return btimes


def ReadStar(file):
	
	f=open(file,'r').readlines()
	
	time=np.empty(len(f))
	flux=np.empty(len(f))
	err=np.empty(len(f))
	sky=np.empty(len(f))
	for i in range(0,len(f)):
		time[i]=float(f[i].split()[0])
		flux[i]=float(f[i].split()[1])
		err[i]=float(f[i].split()[2])
		sky[i]=float(f[i].split()[3])
		
	#time,flux=np.loadtxt(file,usecols=[0,1],unpack=True)
		
	return time,flux,err,sky


def OutputFile(starid,section,time,flux,err,sky):

	name='star_%s_%s_wES_f1.lc.txt' % (starid,section)
	if os.path.isfile==True:
		print "Overwritting old %s file" % (file)
		comm='rm -rf %s' % (file)
		os.system(comm)	
	
	f=open(name,'w')
		
	for i in range(0,len(time)):
		line="%.8f   %.3f   %.3f   %.3f\n" % (time[i],flux[i],err[i],sky[i])
		f.write(line)
	
	f.close()
	
	#z=np.concatenate((time,flux)).reshape(2,len(time)).transpose()
	#np.savetxt(name,z,fmt='%.8f    %.3f')

	return 0

# on starbase the times are exactly the same
# but on the macbook there are a few sig figs missing
# hence like this below its much faster than on the macbook
def Get_n(file,btimes):

	time,flux,err,sky=ReadStar(file)
	print "[Read LC] %s done..." % (file)
		
	flagged=np.zeros(len(time),float)	
		
	for j in range(0,len(btimes)):
		w=np.where(time==btimes[j])
		#for k in range(0,len(time)):
		#	if abs(time[k] - btimes[j]) < 0.000001:
		if len(w[0]) > 0:
			flagged[w[0][0]]=1.0
	
	n=np.where(flagged<1.0)
	
	return n,time,flux,err,sky,flagged
	
	
############################################################
######################## MAIN ##############################
############################################################

# toggles
output = 1

# read in bad times
btimes=GetBadTimes()

# get star list
templist=cmd.getoutput('ls star_*').split('\n')

# get image subsection
section="%s_%s" % (templist[0].split('.')[0].split('_')[2],templist[0].split('.')[0].split('_')[3])


for i in range(0,len(templist)):
	# star id string
	starid=templist[i].split('_')[1]
	
	# get indexed array with good times
	n,time,flux,err,sky,flagged=Get_n(templist[i],btimes)
	
	# output file
	if output > 0:
		out=OutputFile(starid,section,time[n],flux[n],err[n],sky[n])
		if out != 0:
			print "Problem outputing filtered file, exiting!"
			sys.exit()
	
	print "[Filter] %d/%d" % ((i+1),len(templist))



	