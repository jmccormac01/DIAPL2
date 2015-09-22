# ----------------------------------------------------------------------------------
#
# DIAPL_ListBadFrames.py - 	a program to filter out bad frames based on seeing
#							and altitude of observation
#				
#	Outputs a file with HJD of bad frame. This can then be crossed checked to 
#	remove bad frames from period searches using DIAPL_FilterStarList.py
# 
#	


# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 13/11/11 - 	code writen
# 13/11/11 -	code tested with staralt plot, alts and fwhms working fine
#
#


from numpy import * 
import numpy as np
import os, os.path, sys
from math import *
import commands
import pyfits

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

def GetLSThrs(image):

	f=pyfits.open(image)
	t=f[0].header['ST']
	
	hr=float(t.split(':')[0])
	mn=float(t.split(':')[1])
	sc=float(t.split(':')[2])

	lst=hr+(mn/60.0)+(sc/(60.0*60.0))

	return lst


def GetAlt(image):

# original fortran code from ING's staralt program
#
#C ********************************************************
#C * This calculates elevation from equatorial coords.
#C * No atmospheric refraction correction is applied.
#C *
#C * Inputs : Latitude (fi) in degrees, Altitude in meters,
#C *          Local Sidereal Time (lst) in hours,
#C *          RA (alfa, in hours),
#C *          DEC (delta, in degrees)
#C *          Constants: pi,rad (=pi/180),gra (=1/rad),
#C * Outputs: Alt (degrees)
#C *          ALT = 0 -> Horizon       ALT = 90 -> Zenith
#C ********************************************************
#      subroutine elev(lst,fi,altitude,alfa,delta,alt)
#      implicit none
#      double precision delta,fi,lst,alfa,alt
#      double precision altitude,sinh,pi,rad,gra
#      pi = 4.*atan(1.)
#      rad = pi/180.
#      gra = 180./pi
#      sinh = sin(rad*fi)*sin(rad*delta) +
#     &  cos(rad*delta)*cos(rad*fi)*cos(rad*15.*(lst-alfa))
#     if(abs(sinh).ge.1) then
#	 alt = 90.*sinh/abs(sinh)
#     else
#	 alt = gra*asin(sinh)
#     endif
#     end
#
# La Palma 342.1184 E  +28.7606       2326 m
	
	fi=28.7606
	lst=GetLSThrs(image)
	
	# RA and DEC in hours/deg fro M71
	alfa=19.8961417
	delta=18.7784167

	rad = pi/180.0
	gra = 180.0/pi

	sinh = (sin(rad*fi)*sin(rad*delta)) + (cos(rad*delta)*cos(rad*fi)*cos(rad*15.0*(lst-alfa)))
	
	if abs(sinh) > 1.0:
		alt=90.0*sinh/abs(sinh)
		
	else:
		alt=gra*arcsin(sinh)

	return alt


templist=commands.getoutput('ls M71-2011*.fits').split('\n')

frame_fwhms=np.empty(len(templist))
altlist=np.empty(len(templist))

badlist=[]

for j in range(0,len(templist)):

	command1='sex ' + str(templist[j]) 
	os.system(command1)

	# read in the file
	s=Readcoords('test.cat')

	flen=len(s)

	flagged=np.zeros(flen,float)
	num=np.empty(flen)
	x=np.empty(flen)
	y=np.empty(flen)
	flux_auto=np.empty(flen)
	flux_err=np.empty(flen)
	fwhm=np.empty(flen)
	xpeak=np.empty(flen)
	ypeak=np.empty(flen)
	
	
	for i in range(0,flen):
		num[i]=float(s[i].split()[0])
		x[i]=float(s[i].split()[1])
		y[i]=float(s[i].split()[2])
		flux_auto[i]=float(s[i].split()[3])
		flux_err[i]=float(s[i].split()[4])
		fwhm[i]=float(s[i].split()[5])
		xpeak[i]=int(s[i].split()[6])-1
		ypeak[i]=int(s[i].split()[7])-1
	
	
	hdulist = pyfits.open(templist[j])
	ref_data=hdulist[0].data
	
	# flag True the saturated or edge stars
	for i in range(0,flen):
		if ref_data[ypeak[i],xpeak[i]] > 40000.00 or xpeak[i] < 20.0 or xpeak[i] > 1000.0 or ypeak[i] < 20.0 or ypeak[i] > 1000.0:
			flagged[i]=1.0
	
	index=np.where(flagged < 1.0)
	
	# median FWHM for frame rejection
	frame_fwhms[j]=median(fwhm[index])

	# add alt to the header
	altlist[j]=GetAlt(templist[j])
	
	# remove the catalog for the next image
	os.system('rm -rf test.cat')
	
	print "\nImage [%d/%d] " % (j,len(templist))
	print "Median FWHM: %.2f" % (frame_fwhms[j])
	print "Altitude: %.2f" % (altlist[j])
	
	# check if alt or seeing is bad - if so add to list	
	if frame_fwhms[j] > 6.5 or altlist[j] < 35.0:	
		badlist.append(templist[j])
		print "***FRAME REJECTED***"
		
# make a list of bad images with just their
# HJD values to compare in next step

name='batch1-badtimes_6.5_35.lc.txt'

if os.path.isfile(name) == True:
	print "Overwritting old badtimes file!"
	comm='rm -rf %s' % (name)
	os.system(comm)

# open file for bad times
file=open(name,'w')

for i in range(0,len(badlist)):
	hdulist=pyfits.open(badlist[i])
	hjd=hdulist[0].header['HJD-MID']
	
	line="%s\n" % (hjd)
	file.write(line)

file.close()












	

	
	