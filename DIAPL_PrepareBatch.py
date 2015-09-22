# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_PrepareBatch.py - a program to prepare the batches of images 
#						
# 		This gets all of the images for a given batch and copies them to a folder.	
#		This folder is the given in diapl_setup.par as the image dir
#		
#		This also makes the list for diapl
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 14/10/11 - 	code writen
# 14/10/11 -	code tested
# 28/10/11 - 	batch2 setup added				
#


import commands
import os

templist=commands.getoutput('ls').split('\n')

batchdirs=[]
totalimages=0

for i in range(0,len(templist)):
	if str(templist[i][:4]) == '2011':
		
		#	BATCH1
		#
		#if str(templist[i].split('-')[1]) == '06' and int(templist[i].split('-')[2]) > 25:
		#	path=str(templist[i]) + "/reduced/"
		#	os.chdir(path)
		#	t1=commands.getoutput('ls M71-*.fits').split('\n')
		#	print "%s: %d images" % (templist[i], len(t1))  
		#	totalimages=totalimages+len(t1)
		#	os.system('cp M71-*.fits /data/waspr1/jmcc/M71-2/batch1/')
		#	batchdirs.append(templist[i])
		#	os.chdir('../../')
		#
		#if str(templist[i].split('-')[1]) == '07':
		#	path=str(templist[i]) + "/reduced/"
		#	os.chdir(path)
		#	t1=commands.getoutput('ls M71-*.fits').split('\n')
		#	print "%s: %d images" % (templist[i], len(t1))  
		#	totalimages=totalimages+len(t1)
		#	os.system('cp M71-*.fits /data/waspr1/jmcc/M71-2/batch1/')
		#	batchdirs.append(templist[i])
		#	os.chdir('../../')
		#
		#if str(templist[i].split('-')[1]) == '08' and int(templist[i].split('-')[2]) < 30:	
		#	path=str(templist[i]) + "/reduced/"
		#	os.chdir(path)
		#	t1=commands.getoutput('ls M71-*.fits').split('\n')
		#	print "%s: %d images" % (templist[i], len(t1))  
		#	totalimages=totalimages+len(t1)
		#	os.system('cp M71-*.fits /data/waspr1/jmcc/M71-2/batch1/')
		#	batchdirs.append(templist[i])
		#	os.chdir('../../')
		
		#	BATCH 2
		#
		if str(templist[i].split('-')[1]) == '09' and int(templist[i].split('-')[2]) > 1:	
			path=str(templist[i]) + "/reduced/"
			os.chdir(path)
			t1=commands.getoutput('ls M71-*.fits').split('\n')
			print "%s: %d images" % (templist[i], len(t1))  
			totalimages=totalimages+len(t1)
			os.system('cp M71-*.fits /data/waspr1/jmcc/M71-2/batch2/')
			batchdirs.append(templist[i])
			os.chdir('../../')
		
		if str(templist[i].split('-')[1]) == '10':	
			path=str(templist[i]) + "/reduced/"
			os.chdir(path)
			t1=commands.getoutput('ls M71-*.fits').split('\n')
			print "%s: %d images" % (templist[i], len(t1))  
			totalimages=totalimages+len(t1)
			os.system('cp M71-*.fits /data/waspr1/jmcc/M71-2/batch2/')
			batchdirs.append(templist[i])
			os.chdir('../../')
		
print "\nTotal nights: %d" % (len(batchdirs))
print "Total images: %d" % (totalimages)			
			
os.chdir('/data/waspr1/jmcc/M71-2/batch2/')	

imlist=commands.getoutput('ls M71-*.fits').split('\n')
		
name='M71-IMAGES-BATCH2.txt'

FILE=open(name,'w')

for i in range(0,len(imlist)):
	line="%s\n" % (imlist[i])
	FILE.write(line)
	
FILE.close()

print "%d Images written to M71-IMAGES-BATCH2.txt" % (len(imlist))	
	
	
	
		
			