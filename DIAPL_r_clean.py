
# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_r_clean.py - 	a python wrapper for Francesca's r_clean
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
# 02/12/11 -	added short period searching with r_clean2 - SX Phe objects	
# 19/12/11 -	added CLEAN period searching with r_clean3
# 19/05/12 - 	added scargle + sig searching with r_clean4
# 24/05/12 - 	added correct name to fix scargle.d instead of scargle.dat
# 


import commands as cmd
import os, os.path, sys, time
import numpy as np
#from pylab import *

def FixNames():

	t=cmd.getoutput('ls *.d').split('\n')
	
	for i in range(0,len(t)):
		n2="%s.dat" % (t[i].split('.')[0])
		os.rename(t[i],n2)
		print "%s --> %s" % (t[i],n2)
	
	return 0

def EditRClean(image,toggle):
	
	if toggle == 0:
		file='r_clean'
	if toggle == 1:
		file='r_clean2'
	if toggle == 2:
		file='r_clean3'
	if toggle == 4:
		file='r_clean4'
		
	f=open(file,'r')
	s=f.readlines()
	f.close()
	
	s[2]='set fileroot=%s\n' % (image)
	
	if toggle == 0:
		s[16]='%s_logfile.dat\n' % (image.split('.')[0])
		s[35]='%s_scargle.dat\n' % (image.split('.')[0])
	
	if toggle == 1:
		s[16]='%s_logfile.dat\n' % (image.split('.')[0])
		s[37]='%s_scargle.dat\n' % (image.split('.')[0])
	
	if toggle == 2:
		s[17]='%s_logfile.dat\n' % (image.split('.')[0])
		s[38]='%s_clean.dat\n' % (image.split('.')[0])
	
	if toggle == 4:
		s[17]='%s_logfile.dat\n' % (image.split('.')[0])
		s[39]='%s_scargle.dat\n' % (image.split('.')[0])
	
	
	f=open(file,'w')
	
	for i in range(0,len(s)):
		f.write(str(s[i]))
	f.close()
	
	return 0

	
def RunPeriod(toggle):	
	
	t=cmd.getoutput('ls star_0*').split('\n')
	
	for i in range(0,len(t)):
		
		edit=EditRClean(t[i],toggle)
		if edit != 0:
			print "Problem editing r_clean, exiting!"
			sys.exit()
		
		if toggle == 0:
			os.system('./r_clean')
			print "[r_clean] %d/%d..." % ((i+1),len(t))
		if toggle == 1:
			os.system('./r_clean2')
			print "[r_clean2] %d/%d..." % ((i+1),len(t))
		if toggle == 3:
			os.system('./r_clean3')
			print "[r_clean3] %d/%d..." % ((i+1),len(t))
		if toggle == 4:
			os.system('./r_clean4')
			print "[r_clean4] %d/%d..." % ((i+1),len(t))
			
		time.sleep(3)	
				
	return 0


def GetSig():

	t=cmd.getoutput('ls *_logfile.dat').split('\n')
	
	starid,sig=[],[]
	
	for i in range(0,len(t)):
		s=open(t[i]).readlines()
		starid.append(t[i].split('_')[1])
		sig.append(s[8].split()[-1])		

		print "%s Sig1000: %s" % (starid[i],sig[i])

	return 0


def GetCFD(powers):

	powers.sort()
	
	num=np.empty(len(powers))
	
	for i in range(0,len(powers)):
		n=0
		for j in range(0,len(powers)):
			if powers[j]>powers[i]:
				n=n+1
				
		# get n, normalized to the length of the sample
		num[i]=float(n)/len(powers)
	
	return num


def GetPowers():

	t=cmd.getoutput('ls star_*_logfile.dat').split('\n')
	
	section='%s_%s' % (t[0].split('_')[2],t[0].split('_')[3])
	
	powers=np.empty(len(t))
	for i in range(0,len(t)):
		f=open(t[i],'r').readlines()
		powers[i]=float(f[6].split()[-1])

	num=GetCFD(powers)

	# the power corresponding to num=0.01 is the power above which
	# any detection can be confirmed as real above the noise. 
	check=abs(num-0.01)
	loc=np.where(check==min(check))
	fap_power=powers[loc]
	
	print "FAP[0.01]: %f" % (fap_power)
	print "Light curves: %d" % (len(powers))

	# False Alarm Prob - (FAP) - lines
	fap_x=[6,8,10,12,14,16,18,20,22,]
	fap_y=[0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]

	figure(1)
	plot(powers,num,'k-')
	plot(fap_x,fap_y,'r--')
	title("Probability Distribution Function")
	ylabel("Probability of Real Signal")
	xlabel("Periodogram Power")
	xlim(8.0,20.0)
	show()

	return section,powers,num,fap_x,fap_y


def OutputFile(section,powers,num,fap_x,fap_y):
	
	from datetime import date
	
	date=date.today()
	
	name='PDF_shuf_%s_%s.lc.txt' % (section,date)
	name2='PDF_FAP_%s_%s.lc.txt' % (section,date)
	
	if os.path.exists(name) == True:
		print "Overwritting old PDF file!"
		comm="rm -rf %s" % (name)
		os.system(comm)
	if os.path.exists(name2) == True:
		print "Overwritting old PDF FAP file!"
		comm2="rm -rf %s" % (name2)
		os.system(comm2)
	
	z=np.concatenate((powers,num)).reshape(2,len(powers)).transpose()
	np.savetxt(name,z,fmt='%.7f    %.8f')

	z2=np.concatenate((fap_x,fap_y)).reshape(2,len(fap_x)).transpose()
	np.savetxt(name2,z2,fmt='%d    %.2f')

	return 0

# combine the results from the subsections
def Combine(combine_out):
	
	os.chdir('/Volumes/DATA/NITES/Results/M71/Photometry/Batch1/A8.0D12S8C12/')
	section=cmd.getoutput('ls').split('\n')
	
	p=[]
	for i in range(0,len(section)):
		newdir='%s/shuf/' % (section[i])
		os.chdir(newdir)
	
		t=cmd.getoutput('ls PDF_shuf*')
	
		x=np.loadtxt(t,usecols=[0])
		
		for j in range(0,len(x)):
			p.append(x[j])
			
		os.chdir('../../')
	
	num_t=GetCFD(p)
	
	# change powers to an np.array
	powers=np.copy(p)
	
	# the power corresponding to num=0.01 is the power above which
	# any detection can be confirmed as real above the noise. 
	check=abs(num_t-0.01)
	loc=np.where(check==min(check))
	fap_power=powers[loc]
	
	print "FAP[0.01]: %f" % (fap_power)
	print "Light curves: %d" % (len(powers))

	# False Alarm Prob - (FAP) - lines
	fap_x=[6,8,10,12,14,16,18,20,22,]
	fap_y=[0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01]

	figure(1)
	plot(powers,num_t,'k-')
	plot(fap_x,fap_y,'r--')
	title("Probability Distribution Function")
	ylabel("Probability of Real Signal")
	xlabel("Periodogram Power")
	xlim(8.0,20.0)
	show()

	if combine_out > 0:
		section="total"
		out=OutputFile(section,powers,num_t,fap_x,fap_y)
		if out != 0:
			print "Problem outputing PDF file, exiting!"
			sys.exit()	
		
	return num_t

# Main #

# toggle subroutines off/on = 0/1
cor_name=0
r_period=0
get_sig=1
get_powers=0
output=0
combine=0
combine_out=0

# type of period search
toggle=4

if cor_name>0:
	cor=FixNames()
	
	if cor != 0:
		print "Problem correcting names, exiting!"
		sys.exit()

if r_period>0:	
	ran=RunPeriod(toggle)
	
	if ran != 0:
		print "Problem running PERIOD, exiting!"
		sys.exit()

if get_sig>0:
	sig=GetSig()
	
	if sig != 0:
		print "Problem getting Sig1000, exiting!"
		sys.exit()
	
if get_powers>0:
	section,powers,num,fap_x,fap_y=GetPowers()
	
if output>0:
	out=OutputFile(section,powers,num,fap_x,fap_y)
	if out != 0:
		print "Problem outputing PDF file, exiting!"
		sys.exit()

if combine > 0:
	num_t=Combine(combine_out)
	


