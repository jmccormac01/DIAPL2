# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_CheckOutputs.py - 	a program to check the outputs from DIAPL 
#						
# 		Give the number of frames rejected per subsection + a total				
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 18/10/11 - 	code writen
# 18/10/11 -	code tested
#
#	

import commands
import os

def GetNumImagesOriginal():
	
	s=open('M71-IMAGES-BATCH1.txt','r').readlines()
	n=len(s)
	
	return n
	
def GetNumSubtractedImages():

	t=[]
	x=os.getcwd()
	y=os.listdir(x)
	for i in range(0,len(y)):
		if y[i].startswith('s_M71'):
			t.append(y[i])
	
	n=len(t)

	return n
	
def GetNumRImages():
	
	t=commands.getoutput('ls rimages*').split('\n')
	n=[]
	
	for i in range(0,len(t)):
		s=open(t[i],'r').readlines()
		n.append(len(s))
	
	
	return t,n


orig=GetNumImagesOriginal()
origsubs=orig*4
print "\nThere were %d original images (%d subsections)" % (orig, origsubs )

sub=GetNumSubtractedImages()
print "There are %d subtracted image subsections" % (sub)

rfiles,rnum=GetNumRImages()
for i in range(0,4):
	print "In %s there are %d images" % (rfiles[i],rnum[i])
	
q=sum(rnum)
rejected=origsubs-q
frac_rej=(float(rejected)/float(origsubs))*100

print "%d subsections have been rejected! [%f percent]" % (rejected,frac_rej)

