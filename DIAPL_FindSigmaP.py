
#
# Bootstrap to find Sigma_P for M71 survey
# Run PDM on range -0.1 < P < +0.1 d for best fitting period to find error
#

# test on a sine curve to find known period

import numpy as np
from pylab import *
import random as rnd
import os
import commands as cmd
import time

niterations=500

def GetPDMInputs():

	#f=open('/Volumes/DATA/nites/Results/M71/PeriodErrors/BootStrap/BootStrapPeriodList_new.txt').readlines()
	f=open('BootStrapPeriodList_new.txt').readlines()
	
	name=[]
	period=np.empty(len(f))
	for i in range(0,len(f)):
		name.append(f[i].split()[0])
		period[i]=float(f[i].split()[1])
	
	print "Got PDM inputs..."
	
	return name,period

def MakeBootStrapLCs(name):

	for k in range(0,len(name)):
		
		folder=name[k].split('_')[0]
		
		t=os.listdir('.')
		if folder in t:
			print "%s lc's are done..." % (folder)
		
		if folder not in t:
			os.system('mkdir %s' % folder)
		
			os.chdir(folder)
			os.system('mv ../%s .' % name[k])
			
			# load in the real data
			x,y=np.loadtxt(name[k],usecols=[0,1],unpack=True)
			
			for j in range(0,niterations):
				# create new bootstrapped areas to hold randomly selected data
				x_bs=np.empty(len(x))
				y_bs=np.empty(len(y))
				
				# bootstrap the data
				for i in range(0,len(x)):
					r=rnd.randrange(0,len(x),1)
					x_bs[i]=x[r]
					y_bs[i]=y[r]
					
				# sort the data according to time for normalities sake
				temp=zip(x_bs,y_bs)
				temp.sort()
				x_new,y_new=zip(*temp)
				
				# create a new file per bootstrap
				fname="%s_BootStrap_%04d.lc.txt" % (name[k].split('_')[0],j+1)
				z=np.concatenate((x_new,y_new)).reshape(2,len(x_new)).transpose()
				np.savetxt(fname,z,fmt='%.8f    %.5f')
				print "%s saved..." % (fname)
			
			print "[%d/%d] %s done..." % (k+1,len(name),folder)
			
			os.chdir("../")
			time.sleep(2)
		
	return 0

# Vik's PERIOD Program
# run PDM on the bootstrapped lc's 
# coded in looking for existing files to save time by not duplicating them
def RunPDM(v,p):
	
	vdir="%s" % v.split('_')[0]
	os.chdir(vdir)
	
	t1=cmd.getoutput('ls').split('\n')
	if 'r_clean5' not in t1:
		os.system('cp ~/bin/r_clean5 .')
					
	t2=cmd.getoutput('ls v*_BootStrap_*.lc.txt').split('\n')
				
	p_frac=p/1000.0
	
	f_min=1.0/(p+p_frac)
	f_max=1.0/(p-p_frac)
	f_step=((f_max-f_min)/1000.0)/50.0
	print f_step
	time.sleep(2)
	
	# len(t2)
	for i in range(0,1):
				
		# edit the r_clean file
		f=open('r_clean5','r').readlines()
		f[2]="set fileroot=%s\n" % (t2[i])
		f[16]="%s.dat\n" % (t2[i].split('.')[0])		
		f[24]="%.9f\n" % (f_min)
		f[26]="%.9f\n" % (f_max)
		f[28]="%.9f\n" % (f_step)	
			
		f2=open('r_clean5','w')
		for k in range(0,len(f)):
			f2.write(f[k])
		f2.close()
				
		os.system('./r_clean5')
			
	os.chdir('../')	
				
	return 0
	
# query the PDM outputs and get the final error
# print this to file - possibly the input file from the beginning	
def GetPDMError():
	
	t1=cmd.getoutput('ls *.dat').split('\n')
	t2=cmd.getoutput('ls *BootStrap*.lc.txt').split('\n')
	
	if len(t2) != len(t1):
		print "%d dat files and %d lcs, exiting" % (len(t1), len(t2))
		exit()
	
	pdm_stats=np.empty(len(t1))
	periods=np.empty(len(t1))
	
	for i in range(0,len(t1)):
		f=open(t1[i]).readlines()
		
		pdm_stats[i]=float(f[5].split()[-1])
		periods[i]=float(f[6].split()[-1])
	
	print periods
	
	pdm_err=np.std(periods)
	print "PDM error: %.6f" % (pdm_err)
	
	figure(10)
	subplot(121)
	plot(pdm_stats,'r-')
	subplot(122)
	plot(periods,'g-')
	
	show()
	
	return periods, pdm_stats

# need a copy of the login.cl for IRAF to work
def LoadIRAF():
	# check for IRAF login file
	if os.path.exists('login.cl') == False:
		os.system('cp ~/login.cl .')

	from pyraf import iraf
	x=raw_input('IRAF loaded: Press RETURN')
	
	return iraf

# IRAF
def PDM_IRAF(v,p):

	p_frac=p/1000.0
	
	p_min=p-p_frac
	p_max=p+p_frac
	
	vdir="%s" % v.split('_')[0]
	print "%s" % (vdir)
	os.chdir(vdir)
	
	t3=cmd.getoutput('ls ').split('\n')
	
	make_pdm_files = 1
	
	for i in range(0,len(t3)):
		if t3[i].startswith('pdmbatch') == True:
			print "PDM files exist, using them..."
			make_pdm_files = 0
			break
		
	if make_pdm_files > 0:
		t=cmd.getoutput('ls v*_BootStrap_*.lc.txt').split('\n')
			
		for i in range(0,len(t)):
			n1="pdmmeta_%d" % int(i+1)
			n2="pdmbatch_%d" % int(i+1)
			iraf.pdm(infiles=t[i],metafil=n1, batchfi=n2, interac="no", minp=p_min, maxp=p_max)
			print "[%d/%d]" % (int(i+1),len(t))

	t2=cmd.getoutput('ls pdmbatch*').split('\n')

	periods=np.empty(len(t2))
	
	for j in range(0,len(t2)):
		f=open(t2[j]).readlines()
		periods[j]=float(f[2].split()[2].split(',')[0])
		
	err=np.std(periods)
	
	print "%s err: %.6f" % (vdir, err)
		
	return periods

					
# MAIN

# toggles
make_lcs = 0
run_pdm_vik = 0
get_pdm_err = 0

# use this one!
run_pdm_iraf = 1

name,period=GetPDMInputs()

# make boot strapped lcs
if make_lcs > 0:
	made=MakeBootStrapLCs(name)

# run pdm
if run_pdm_vik > 0:
	for i in range(0,1):
		ran=RunPDM(name[i],period[i])
		
	if get_pdm_err > 0:
		periods,pdm_stats=GetPDMError()

if run_pdm_iraf > 0:
	iraf=LoadIRAF()
	
	iraf.noao(_doprint=0)
	iraf.astutil(_doprint=0)
	
	for i in range(4,5):
		periods=PDM_IRAF(name[i],period[i])
