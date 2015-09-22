# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_MakeSourceList.py - 	a program to make a source list from template for 
#								DIAPL2 using sextractor.
# 						
#		This script prepares the file for the 'phot.bash' script of DIAPL2
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 12/10/11 - 	code writen
# 12/10/11 -	code tested
# 14/10/11 -	added cps to the program	
# 18/10/11 - 	updated regions for M71-batch1
# 28/10/11 - 	added *.coords files for PSF fitting routine to get magnitudes
# 16/07/12 -	Updated file names to be general after M71 survey was finished
# 27/05/14 - 	Made tplr*_* images generic sizes
#

from numpy import * 
import numpy as np
import os, sys
import commands
import pyfits, time

# read ascii file
def Readcoords(ascifile):                              
	f=file(ascifile,'r')
	r1=f.readline()
	r2=f.readline()
	r3=f.readline()
	r4=f.readline()
	r5=f.readline()
	r6=f.readline()
	r7=f.readline()
	r8=f.readline()
	s=f.readlines()

	return s

def Cps():
	
	os.system('cp ~/default.* .')
	
	return 0


copied=Cps()
if str(copied) != '0':
	print "Problem running 'cps', exiting!"
	exit()
	
# run Sextractor on the first frame in the directory.
templist=commands.getoutput('ls tplr*.fits').split('\n')

# get object name
objects=commands.getoutput('ls *m.fits').split('\n')
object=objects[0].split('-')[0]

for i in range(0,len(templist)):
	
	# run sextractor
	command='sex ' + str(templist[i]) 
	os.system(command) 
		
	# read in the file
	s=Readcoords('test.cat')
	
	# data arrays
	flagged=np.zeros(len(s),float)
	num=np.empty(len(s))
	x=np.empty(len(s))
	y=np.empty(len(s))
	flux_auto=np.empty(len(s))
	flux_err=np.empty(len(s))
	fwhm=np.empty(len(s))
	xpeak=np.empty(len(s))
	ypeak=np.empty(len(s))
	
	# split up the sextractor output
	for j in range(0,len(s)):
		if str(s[j][0]) != '#':
			num[j]=int(s[j].split()[0])
			x[j]=float(s[j].split()[1])
			y[j]=float(s[j].split()[2])
			flux_auto[j]=float(s[j].split()[3])
			flux_err[j]=float(s[j].split()[4])
			fwhm[j]=float(s[j].split()[5])
			xpeak[j]=int(s[j].split()[6])-1
			ypeak[j]=int(s[j].split()[7])-1
	
	hdulist = pyfits.open(templist[i])
	ref_data=hdulist[0].data
	
	print "Have you checked the image limits?"
	yn=raw_input("(e.g. y): ")
	if str(yn) != 'y':
		print "Check limits, exiting..."
		sys.exit()
	
	
	# flag '1.0' the saturated or edge stars
	if str(templist[i]) == 'tplr1_1.fits':
		xs,ys=pyfits.open('tplr1_1.fits')[0].data.shape
		for j in range(0,len(s)):
			if ref_data[ypeak[j],xpeak[j]] > 45000.00 or xpeak[j] < 25.0 or ypeak[j] < 25.0 or xpeak[j] > xs-5 or ypeak[j] > ys-5:
				flagged[j]=1.0
				
	if str(templist[i]) == 'tplr1_2.fits':
		xs,ys=pyfits.open('tplr1_2.fits')[0].data.shape
		for j in range(0,len(s)):
			if ref_data[ypeak[j],xpeak[j]] > 45000.00 or xpeak[j] < 25.0 or ypeak[j] < 5.0  or xpeak[j] > xs-5 or ypeak[j] > ys-25:
				flagged[j]=1.0
				
	if str(templist[i]) == 'tplr2_1.fits':
		xs,ys=pyfits.open('tplr2_1.fits')[0].data.shape
		for j in range(0,len(s)):
			if ref_data[ypeak[j],xpeak[j]] > 45000.00 or xpeak[j] < 5.0  or ypeak[j] < 25.0 or xpeak[j] > xs-25 or ypeak[j] > ys-5:
				flagged[j]=1.0
					
	if str(templist[i]) == 'tplr2_2.fits':
		xs,ys=pyfits.open('tplr2_2.fits')[0].data.shape
		for j in range(0,len(s)):
			if ref_data[ypeak[j],xpeak[j]] > 45000.00 or xpeak[j] < 5.0  or ypeak[j] < 5.0  or xpeak[j] > xs-25 or ypeak[j] > ys-25:
				flagged[j]=1.0
		
	
	# index the good stars
	index=np.where(flagged < 1.0)
	
	# compare the original numbers of stars to 
	# those after removing non-linear etc
	print "\nOriginal: " + str(len(s)) + " stars"
	print "Now: " + str(len(xpeak[index])) + " stars"
	
	# make the names of the new files
	name=str(object) + '_' + str(templist[i].split('.')[0].split('r')[1]) + '.coo'
	name2=str(object) + '_' + str(templist[i].split('.')[0].split('r')[1]) + '.flux'
	
	# check the names of the new files
	print str(templist[i]) + ": " + str(name) + ": " + str(name2)
	
	# wait a sec
	time.sleep(1)

	# make the *.coo files		
	FILE=open(name,"w")
		
	for j in range(0,len(num[index])):
		linept2="%.2f   %.2f\n" % (x[index][j],y[index][j])
		line = str(j+1) + "   " + str(linept2)
		print line
		FILE.write(line)
			
	FILE.close		
	
	# make the *.flux files	
	FILE2=open(name2,"w")

	for j in range(0,len(num[index])):
		line2pt2="%.2f   %.2f   %.2f   %.2f\n" % (x[index][j],y[index][j],flux_auto[index][j],flux_err[index][j])
		line2 = str(j+1) + "   " + str(line2pt2)
		FILE2.write(line2)
			
	FILE2.close	 
	
	print "Do you want a tvmark file?"
	tvmark_yn=raw_input("(e.g. y): ")
	if str(tvmark_yn) == 'y':
		
		name3=str(object) + '-tvmark_' + str(templist[i].split('.')[0].split('r')[1]) + '.coo'
		filler=0.5000
		
		FILE3=open(name3,"w")
		
		for i in range(0,len(num[index])):
			line3 = str(x[index][i]) + "   " + str(y[index][i]) + "   " + str(filler) + "   " + str(filler) + "   " + str(filler) + "   " + str(num[index][i]) + "\n" 
			FILE3.write(line3)

		FILE3.close
		
	# remove the cataloge from sextractor
	os.system('rm test.cat')
	