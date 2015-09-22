
# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_TrimImages.py	-	A program to trim images to 1024x1024 for DIAPL
#							
#							*Code to be appended to NITES_pipe.py* 				
#				
# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 27/09/11	- 	code writen
# 27/09/11	-	code tested
#	

import commands 
import os

def Cpl():

	os.system('cp ~/login.cl .')

	return 0
	

# trim images
def TrimImages():
	
	print "\nTrimming images...\n"
	
	templist=commands.getoutput('ls *.fit').split('\n')
	
	for i in range(0,len(templist)):
		image=str(templist[i]+"[2]")
		
		image2="%s_t.fit" % (image.split('.fit')[0])
		
		
		# trimming other images
		#image=str(templist[i]+"[1][258:770,16:528]")
		#image2="%s_t.fits" % (templist[i].split('.')[0])

		iraf.imcopy(input=image,output=image2)
	
		print image + "[" + str(i+1) + "/" + str(len(templist)) + "]"
				
	return 0

f1=Cpl()

from pyraf import iraf
	
x=TrimImages()
if str(x) != "0":
	print "Problem trimming images, exiting!"
	exit()
	
