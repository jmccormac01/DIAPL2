

# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_RepairTruncated.py - 	a python program to repair truncated images
#						
# 						

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 20/12/11 - 	code writen - Happy Birthday!
# 20/12/11 -	code tested
#				
#

import commands as cmd
import pyfits as py
import os, os.path

if os.path.exists('login.cl') == False:
	os.system('cp ~/login.cl .')
	
from pyraf import iraf 

t=cmd.getoutput('ls s_*.fits').split('\n')

for i in range(0,len(t)):
	hdulist=py.open(t[i])
	data=hdulist[0].data
	
	x,y=data.shape
	
	print data.shape
	
	print "Making copy of truncated file..."
	
	image = "%s[1:%d,1:%d]" % (t[i],x,y)
	image2 = "%s-t.fits" % (t[i].split('.')[0])
	
	iraf.imcopy(input=image,output=image2)
	