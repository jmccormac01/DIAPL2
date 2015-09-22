# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_ListImages.py	-	A program to list all the images to be analysed by DIAPL
#							
#							*Code to be appended to NITES_pipe.py* 				
#				
# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 13/09/11	- 	code writen
# 13/09/11	-	code tested
# 22/06/12 	- 	Code changed to ask for date input instead of trying to find it	
#

import commands
import os

def MakeImageList():

	# get list of file names
	templist=commands.getoutput('ls *-*-*.fits').split('\n')
	
	# create a txt file with one file per line
	date=raw_input('Observation Date (e.g. 20120620): ')
	name="%s-%s-Images.txt" % (templist[0].split('-')[0],date)
		
	FILE=open(name,'w')
	
	for i in range(0,len(templist)):
		line = str(templist[i]) + "\n"
		FILE.write(line)
			
	FILE.close()
	
	return 0

x=MakeImageList()

