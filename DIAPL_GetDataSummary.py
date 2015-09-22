
# ----------------------------------------------------------------------------------
#
# DIAPL_GetDataSummarry.py - a program to get a summary of all the data for thesis
#						
#
#	


# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 14/11/11 - 	code writen
# 14/11/11 -	code tested
#
#

import commands
import os
import pyfits

templist=commands.getoutput('ls').split('\n')

im_total = 0
hours_total = 0


for i in range(0,len(templist)):
	if templist[i][:4] == '2011':
		if templist[i][5:7] == "06" or templist[i][5:7] == "07" or templist[i][5:7] == "08":
			file="%s/reduced" % (templist[i])
			os.chdir(file)
			
			imlist=commands.getoutput('ls M71-*.fits').split('\n')
			
			hdulist=pyfits.open(imlist[0])
			exp=hdulist[0].header['EXPTIME']
			min_am=hdulist[0].header['AIRMASS']
			
			hdulist2=pyfits.open(imlist[-1])
			max_am=hdulist2[0].header['AIRMASS']
			
			print "%s - Images: %d Exp: %d Airmass: %s - %s" % (templist[i], len(imlist), int(exp), min_am, max_am)
			
			im_total = im_total + len(imlist)
			hours_total = hours_total + ((len(imlist)*exp)/(60*60))
			
			os.chdir('../../')
			
print "Totals:"
print "\tImages: %d" % (im_total)
print "\tExptime: %.2f" % (hours_total) 