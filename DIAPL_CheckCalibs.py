# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_CheckCalibs.py - a program to check the calibration frames
#						
# 		This gets all of the calibs files if needed and then displays them
#		for each night.				
#
#		Run this script in the top level directory 
#		/Volumes/DATA/NITES/Results/M71/CheckCalibs/
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 14/10/11 - 	code writen
# 14/10/11 -	code tested
#				
#				
	

import commands
import os, os.path
import sys


def GetFiles():
	
	os.chdir('/Volumes/DATA/NITES/Results/M71/CheckCalibs/')
	
	templist=commands.getoutput('ssh jmcc@starbase.mp.qub.ac.uk ls /data/waspr1/jmcc/M71-2/').split('\n')
	
	for i in range(0,len(templist)):
		if str(templist[i][:4]) == '2011':
			if os.path.exists(templist[i]) == False:
				os.mkdir(templist[i])
				os.chdir(templist[i])
			if os.path.exists(templist[i]) == True:
				os.chdir(templist[i])
			
			# flat
			command1='scp -r jmcc@starbase.mp.qub.ac.uk:/data/waspr1/jmcc/M71-2/%s/reduced/calibs/Flat.fits .' % (templist[i])
			os.system(command1)
			
			# bias
			command2='scp -r jmcc@starbase.mp.qub.ac.uk:/data/waspr1/jmcc/M71-2/%s/reduced/calibs/Zero.fits .' % (templist[i])
			os.system(command2)
			
			# dark
			command3='scp -r jmcc@starbase.mp.qub.ac.uk:/data/waspr1/jmcc/M71-2/%s/reduced/calibs/Dark.fits .' % (templist[i])
			os.system(command3)
			
			os.chdir('../')
			
	return 0
	
def CheckFiles():

	os.chdir('/Volumes/DATA/NITES/Results/M71/CheckCalibs/')
	os.system('cp ~/login.cl .')
	
	from pyraf import iraf
	
	templist=commands.getoutput('ls').split('\n')

	flat,dark,bias,night=[],[],[],[]

	for i in range(0,len(templist)):
		os.chdir(templist[i])
		
		print "\nCurrent Directory: " + str(templist[i])
		 
		iraf.display(image='Flat.fits', frame=1)
		#iraf.imexam(frame=1,image='Flat.fits',input='Flat.fits')
		x=raw_input("(Good? y/n): ")
		flat.append(x)
		
		iraf.display(image='Dark.fits', frame=2)
		#iraf.imexam(frame=2,image='Dark.fits',input='Dark.fits')
		y=raw_input("(Good? y/n): ")
		dark.append(y)
		
		iraf.display(image='Zero.fits', frame=3)
		#iraf.imexam(frame=3,image='Zero.fits',input='Zero.fits')
		z=raw_input("(Good? y/n): ")
		bias.append(z)
		
		night.append(templist[i])
		
		os.chdir('../')
	
	return 0
	
	
# main
print "Get calibtration files from starbase?"
yn=raw_input("(e.g. y): ")
if str(yn)=='y':
	f1=GetFiles()
	if f1 != 0:
		print "Problem getting calibs from starbase, exiting!"
		sys.exit()	
	
print "Check the calibs files for each night?"
yn2=raw_input("(e.g. y): ")
if str(yn2)=='y':
	f2=CheckFiles()
	if f2 != 0:
		print "Problem checking calibs, exiting!"
		sys.exit()








