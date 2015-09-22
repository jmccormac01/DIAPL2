# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_SplitPhot.py - 	a program to split the phot output of diapl2 
#						
# 						
#		This script prepares the phot output like that from multphot in IRAF
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 12/10/11	- 	code writen
# 12/10/11	-	code tested
# 24/10/11	- 	added print statements to follow progress
# 24/10/11	-	had to restructure the program due to HUGE *.db files
# 20/07/12	- 	changed M71 dependence
#
#

import os, commands, sys
import pyfits
import fileinput
from numpy import *

def TidyFiles():

	os.mkdir('1_1')
	os.mkdir('1_2')
	os.mkdir('2_1')
	os.mkdir('2_2')

	os.mkdir('original')

	os.system('mv *1_1.txt 1_1')
	os.system('mv *1_2.txt 1_2')
	os.system('mv *2_1.txt 2_1')
	os.system('mv *2_2.txt 2_2')
	
	os.system('mv *.db* original')
	
	print "Files tidied, done!"
	
	return 0


templist=commands.getoutput('ls *_*.db').split('\n')

for q in range(0,len(templist)):

	imlist=[]
	lenfile=0
	i=0
	sub=str(templist[q].split('.')[0].split('_')[1]) + "_" + str(templist[q].split('.')[0].split('_')[2])
	
	for line in fileinput.input([templist[q]]):
		if str(line[0]) == '0':
			
			if i>1:
				FILE.close()
				print "%s written..." % (name)
				
			i=i+1
			imlist.append(line.split()[1])	
			name = 'coords%05d_%s.txt' % (i,sub)
			FILE=open(name,'w')
		
		if str(line[0]) != '0':
			FILE.write(line)
			
		
		lenfile=lenfile+1
	
	# no. stars (real number is -1 but this works better in the line below)
	nstars=lenfile/len(imlist)
	
	print "\nThere are " + str(len(imlist)) + " images, with " + str(nstars-1) + " stars each"
	
	name2="%s-imlist-%s.txt" % (templist[0].split('_')[0],sub)
	FILE2=open(name2,"w")
	for d in range(0,len(imlist)):
		line2=str(imlist[d]) + "\n"
		FILE2.write(line2)
	FILE2.close()	
	
	print "\n%s written..." % (name2)	
		
tidied=TidyFiles()
if tidied != 0:
	print "Problem tidying files up, exiting!"
	sys.exit()
	
	
		
		