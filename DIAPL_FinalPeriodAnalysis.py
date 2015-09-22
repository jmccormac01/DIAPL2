
# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_FinalPeriodAnalysis.py - 	a python interpreter for PERIOD:* results
#						
# 						
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 28/11/11 - 	code writen
# 28/11/11 -	code tested
#				
#

import commands as cmd
import numpy as np
from datetime import date
import os,os.path
import sys

def FixNames():

	t=cmd.getoutput('ls *_sca').split('\n')
	
	for i in range(0,len(t)):
		new_name="%srgle.dat" % (t[i])
		os.rename(t[i],new_name)

	return 0

# toggles
fix_names = 0
outfiles = 0
clean = 1
scargle =0

if fix_names > 0:
	fix=FixNames()
	if fix != 0:
		print "Problem fixing names, exiting!"
		sys.exit()

t=cmd.getoutput('ls *_logfile.dat').split('\n')
t2=cmd.getoutput('ls *_clean.dat').split('\n')

# arrays for peak powers
peak_power=np.empty(len(t))
peak_period=np.empty(len(t))
peak_err=np.empty(len(t))
starid=np.empty(len(t))

# fill peak power arrays
for i in range(0,len(t)):
	starid[i]=t[i].split('_')[1]
	
	f=open(t[i],'r').readlines()
	
	if scargle > 0:
		peak_power[i]=f[6].split()[-1]
		peak_period[i]=f[7].split()[-1]
		peak_err[i]=f[8].split()[-1]
	
	if clean > 0:
		peak_power[i]=f[5].split()[-1]
		peak_period[i]=f[6].split()[-1]
		peak_err[i]=f[7].split()[-1]
	
	print "[Get Peaks] %d/%d" % ((i+1),len(t)) 

if scargle > 0:	
	
	# arrays for Sg
	Sg=np.empty(len(t2))
	real=np.zeros(len(t2))
	
	# get Sg
	for i in range(0,len(t2)):
		freqs,powers=np.loadtxt(t2[i], usecols=[0,1], unpack=True)
		
		Sg[i]=(peak_power[i]-np.average(powers))/(np.std(powers))
		
		if peak_power[i] > (16*Sg[i]):
			 real[i]=1.0
		if peak_power[i] > (10*Sg[i]) and peak_power[i] <= (16*Sg[i]):
			real[i]=2.0
	
		print "[Get Sg] %d/%d" % ((i+1),len(t2)) 	
	
	# print summary	
	n=np.where(real==1.0)
	print "Stars with Power > 16Sg:\n"
	for i in range(0,len(n[0])):
		print "Star %d: P = %f [%f]" % (starid[n[0][i]],peak_period[n[0][i]],peak_power[n[0][i]])
	
	n2=np.where(real==2.0)
	print "Stars with power > 10Sg:\n"
	for i in range(0,len(n2[0])):
		print "Star %d: P = %f [%f]" % (starid[n2[0][i]],peak_period[n2[0][i]],peak_power[n2[0][i]])


if outfiles > 0:
	# make 2 output files
	# 1 for significant detections
	# 1 with the significance of all detections and its factor (i.e. how many Sg)
	section = "%s_%s" % (t[0].split('_')[2],t[0].split('_')[3])
	name = "variables_%s.txt" % (section)
	name2 = "significance_%s.txt" % (section)
	
	if len(n[0]) > 0:
		if os.path.exists(name) == True:
			print "Overwritting %s!" % (name)
			comm="rm -rf %s" % (name)
			os.system(comm)
	
		file=open(name,'w')
		for i in range(0,len(n[0])):
			line = "Star %d: P = %f [%f]\n" % (starid[n[0][i]],peak_period[n[0][i]],peak_power[n[0][i]])
			file.write(line)
	
		file.close()
	
	if os.path.exists(name2) == True:
		print "Overwritting %s!" % (name2)
		comm2="rm -rf %s" % (name2)
		os.system(comm2)
	
	file2=open(name2,'w')
	for i in range(0,len(t2)):
		line2 = "%d   %.6f   %.6f   %.6f   %.3f\n" % (starid[i],peak_period[i],peak_power[i],Sg[i],(peak_power[i]/Sg[i]))
		file2.write(line2)
	
	file2.close()
	
