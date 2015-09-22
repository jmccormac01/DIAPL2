# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_MakeFitsTable.py - 	a program to get RA and DEC from X & Y coords
#						
# 		*Code added to DIAPL_PlotStars.py*				
#
# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 17/10/11 - 	code writen
# 17/10/11 -	code tested with Gaia, getting RA and DEC with 0.5 arcsec
#		
#	

import pyfits 
from pyfits import *  
import numpy as np
from pylab import *
import sys, os, os.path

def GetXY():

	f=open('Coords0.txt').readlines()
	
	x,y=[],[]
	
	for i in range(0,len(f)):
		x.append("%.3f" % (float(f[i].split()[0])))
		y.append("%.3f" % (float(f[i].split()[1])))

	return x,y


def MakeFitsTable(x,y):
	c1=Column(name='x', format='E', array=x)
	c2=Column(name='y', format='E', array=y)
	
	tbhdu=new_table([c1,c2])
	print tbhdu.header.ascardlist() 

	name ='StarPosXY.fits'
	
	if os.path.isfile(name) == True:
		print "\nOverwritting old StarPoxXY file..."
		os.system('rm -rf StarPosXY.fits')
		
	tbhdu.writeto(name)
	
	return 0


def WCS_xy2rd():
	
	if os.path.isfile('StarPosRaDec.fits') == True:
		print "\nOverwritting old StarPoxRaDec file..."
		os.system('rm -rf StarPosRaDec.fits')
	
	os.system('wcs-xy2rd -w wcs.fits -i StarPosXY.fits -o StarPosRaDec.fits')

	return 0


def GetRaDec():
	
	RA,DEC=[],[]
	
	t=pyfits.open('StarPosRaDec.fits')
	tbdata=t[1].data
	
	for i in range(0,len(tbdata)):
		ra=tbdata[i][0]
		dec=tbdata[i][1]
	
		ra1=(ra/15)
		ra2=(fmod(ra1,1)*60)
		ra3=(fmod(ra2,1)*60)
		
		if len(str(ra3).split('.')[0]) < 2:
			ra3="0"+str(ra3)
		
		ratot="%02d:%02d:%s" % (int(ra1),int(ra2),str(ra3)[:5])
		RA.append(ratot)
		
		
		dec1=dec
		dec2=(fmod(dec1,1)*60)
		dec3=(fmod(dec2,1)*60)

		if len(str(dec3).split('.')[0]) < 2:
			dec3="0"+str(dec3)
		
		dectot="%02d:%02d:%s" % (dec1,dec2,str(dec3)[:5])
		DEC.append(dectot)
		
		
	return RA,DEC
	

x,y=GetXY()
tab=MakeFitsTable(x,y)
if tab != 0:
	print "Problem making fits table, exiting!"
	sys.exit()

conv=WCS_xy2rd()
if conv != 0:
	print "Problem converting data, exiting!"
	sys.exit()

RA,DEC=GetRaDec()

for i in range(0,len(x)):
	
	print "X: %.2f Y: %.2f - RA: %s DEC: %s" % (float(x[i]), float(y[i]), RA[i], DEC[i])


