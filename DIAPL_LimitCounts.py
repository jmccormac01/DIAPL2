# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_LimitCounts.py	-	A program to mask negative values in reduced frames
#							
#							*Code to be appended to NITES_pipe.py* 				
#				
# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 13/09/11	- 	code writen
# 13/09/11	-	code tested
#	

# remember pyfits uses [y1:y2,x1:x2] - 0 indexed

import commands
import pyfits
import time

def FixPixValues():

	# get list of file names
	templist=commands.getoutput('ls *-*-*.fits').split('\n')
	
	# edit the values < 1 or > 65000  
	for i in range(0,len(templist)):
		hdulist=pyfits.open(templist[i],mode='update')
		prihdr=hdulist[0].header
		scidata=hdulist[0].data
		
		print "Checking pixels in " + str(templist[i])
		
		y=scidata.shape[0]
		x=scidata.shape[1]
		
		n_pix = y*x
		neg=0
		sat=0
		for j in range(0,y):
			for k in range(0,x):
				if scidata[j,k] > 65000.0:
					scidata[j,k] = 65000.0
					sat=sat+1
				if scidata[j,k] < 1:
					scidata[j,k] = 5
					neg=neg+1
		
		t_now=time.ctime()
		line='Fixed pixels on '+ str(t_now)
		prihdr.add_history(line)
		hdulist.close()
		
		print "\tNegative:\t" + str(neg) + "/" + str(n_pix) + " [" + str((neg/n_pix)*100) + "%]"
		print "\tSaturated:\t" + str(sat) + "/" + str(n_pix) + " [" + str((sat/n_pix)*100) + "%]"
		print  "%s: Pixels Fixed [%d/%d]" % (templist[i],i+1,len(templist))
					
	return 0
	
x = FixPixValues()
	
	