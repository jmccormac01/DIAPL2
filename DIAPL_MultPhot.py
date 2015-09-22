
# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_MultPhot.py - 	a program to extract the photometry for multiple apertures
#						
# 						
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 14/10/11 - 	code writen
# 14/10/11 -	code tested
#				added to copy phot.par to phot dirs to make sure what aperture 
#				was used for the run
# 16/07/12 - 	tidied up a little after M71 survey finished
#
#	

import os
import os.path
import sys

def EditPhotPar(num):
	
	f=open('phot.par','r')
	s=f.readlines()
	f.close()
	
	s[0]='APRAD     = %.1f         // Aperture radius\n' % (num)
		
	f=open('phot.par','w')
	
	for i in range(0,len(s)):
		f.write(str(s[i]))
	f.close()
	
	return 0


def MakeDir(num):
	
	# make new dir in the current workign directory
	newdir="A%.1fD12S8C12" % (num)
	
	if os.path.exists(newdir) == False:
		os.mkdir(newdir)
	
	return 0, newdir


def MoveFiles(ndir):
	
	command1='mv photimages* ' + str(ndir) + '/'
	command2='mv *.db* ' + str(ndir) + '/'
	command3='mv phot_err.log ' + str(ndir) + '/'
	command4='cp phot.par ' + str(ndir) + '/'
	
	x1=os.system(command1)
	x2=os.system(command2)
	x3=os.system(command3)
	x4=os.system(command4)
	
	return 0


def MakeNewPhotLog():
	
	os.system('touch phot_err.log')
	
	return 0
	
		
# Main

# aperture list
aps=[6.0,8.0,10.0]

for i in range(0,len(aps)):
	
	# edit phot file
	f1=EditPhotPar(aps[i])
	if f1 != 0:
		print "Problem editing 'phot.par', exiting!"
		sys.exit()
	
	# run the phot.bash script
	f2=os.system('./phot.bash > phot_err.log')
	if f2 != 0:
		print "Problem running 'phot.bash' script, exiting!"
		sys.exit()

	# make phot dir
	f3,ndir=MakeDir(aps[i])
	if f3 != 0:
		print "Problem making directory for aps[%.1f], exiting!" % (aps[i])
		sys.exit()

	# move files
	f4=MoveFiles(ndir)
	if f4 != 0:
		print "Problem moving files to aps[%.1f] directory, exiting!" % (aps[i])
		sys.exit()

	f5=MakeNewPhotLog()
	if f5 != 0:
		print "Problem making new 'phot_err.log' file, exiting!" 
		sys.exit()
		
	