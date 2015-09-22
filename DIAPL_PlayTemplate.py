# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_PlayTemplate.py	-	A program to display the images to be used for 
#							the template image
#										
#							Used to check for bad frames in the template


# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 10/10/11	- 	code writen
# 10/10/11	-	code tested
# 10/10/11	- 	now takes old list, gets images, displays them, askes for input
#				on them and makes a list of the ones to keep	
# 11/10/11 	-	Added send/get template.list to/from mars	
# 15/07/12	- 	Added fix to importing IRAF problem on LION	
#				*IRAF must be called with a raw_input after to let PYTHON load,
#				then iraf must be returned from the fuction or it is not known
#				about*
#
#	TODO
#				Remove the images from the current directory
#

import time
import os, os.path
import sys

def LoadIRAF():
	# check for IRAF login file
	if os.path.exists('login.cl') == False:
		os.system('cp ~/login.cl .')

	from pyraf import iraf
	x=raw_input('IRAF loaded: Press RETURN')
	
	return iraf

def Readcoords(ascifile):
	f=file(ascifile,'r')
	s=f.readlines()   
		
	return s

def GetTemplateList():
	os.system('scp -r jmcc@mars.ing.iac.es:/data/jmcc/diapl/template.list .')
	
	return 0

def SendTemplateList():
	os.system('scp -r template.list jmcc@mars.ing.iac.es:/data/jmcc/diapl/')
	
	return 0

# load IRAF
iraf=LoadIRAF()

# get the template list file from diapl2 working directory on mars
if os.path.exists('template.list') == False:
	gtfile=GetTemplateList()
	if gtfile != 0:
		print "Problem getting template.list, exiting!"
		sys.exit()


# get the list of template images
s=Readcoords('template.list')

slist=[]
for i in range(0,len(s)):
	slist.append(s[i].split()[0])
	
print "\nThere are " + str(len(slist)) + " images\n"

# get the images from mars
for i in range(0,len(slist)):
		
	if os.path.isfile(slist[i]) == True:
		print str(slist[i]) + " exists, skipping..."
	elif os.path.isfile(slist[i]) == False: 
		command="scp -r jmcc@mars.ing.iac.es:/data/jmcc/diapl/"+str(slist[i]) + " ."
		os.system(command)
		
	
yn=[]
print "\nSelect images to keep or exclude!\n"

# display the images
for i in range(0,len(slist)):
	print "Image: " + str(slist[i]) + " [" + str(i+1) + "/" + str(len(slist)) + "]"  
	iraf.display(image=slist[i],frame='1') 
	q=raw_input("(e.g. y): ")
	yn.append(q)


# enter y to keep a frame, anything else to reject it
new_tpl_list=[]
for i in range(0,len(yn)):
	if str(yn[i]) == 'y':
		new_tpl_list.append(slist[i])


# rename the template list
print "Make new template.list?"
tpl_yn=raw_input("(e.g. y): ")
if str(tpl_yn) == 'y':		
	
	os.system('mv template.list old_template.list')
	
	# write the new template.list
	name = 'template.list'
	
	FILE=open(name,'w')
	
	for i in range(0,len(new_tpl_list)):		
		line=str(new_tpl_list[i]) + "\n"
		FILE.write(line)
	
	FILE.close


# send it?
print "Send template.list to mars?"
send_yn=raw_input("(e.g. y): ")
if str(send_yn) == 'y':
	# send the updated template.list file back to mars
	sntfile=SendTemplateList()
	if sntfile != 0:
		print "Problem sending template.list, exiting!"
		exit()




