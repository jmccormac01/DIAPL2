# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_BLS.py - 	a python wrapper for initbls.pro
#
#		Runs BLS fit to data from DIAPL2 (Use this not BLS.py)
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 31/12/11 - 	code writen - Happy New Year!
# 31/12/11 -	code tested
#				
#

import commands as cmd
import numpy as np
import time, os, os.path
import pylab as pl

def EditBLSScript(star):
	
	f=open('initbls_py.pro').readlines()
	
	f[8]='readcol,"%s",t,f\n' % (star)
	f[23]='sid=%s\n' % (star.split('.')[2].split('_')[1])
	f[33]='openw, lun, "%s.%s.%s_bls.dat", /get_lun\n' % (star.split('.')[0], star.split('.')[1],star.split('.')[2])
	
	f2=open('initbls_py.pro','w')
	for i in range(0,len(f)):
		f2.write(f[i])
	f2.close()
	
	return 0


def GetStartsEnds(time):
	
	starts=[0]
	ends=[]
	
	for j in range(0,len(time)):
		if j < len(time)-1:
			if time[j+1] - time[j] > 0.5:
				starts.append(j+1)
				ends.append(j)
	
	ends.append(len(time))
	
	return starts,ends


def BinLcs(t):

	for i in range(0,len(t)):
		time,flux=np.loadtxt(t[i],usecols=[0,1],unpack=True)
		starts,ends=GetStartsEnds(time)
	
		time_binned=np.array([])
		flux_binned=np.array([])
		
		# bin in ~5min bins x10
		for j in range(0,len(starts)):
			time_b=np.empty((ends[j]-starts[j])/10)
			flux_b=np.empty((ends[j]-starts[j])/10)
			for k in range(0,(ends[j]-starts[j])/10):
				time_b[k]=np.average(time[starts[j]:ends[j]][(k*10):(10*(k+1))])
				flux_b[k]=np.average(flux[starts[j]:ends[j]][(k*10):(10*(k+1))])
				
			time_binned=np.hstack((time_binned,time_b))
			flux_binned=np.hstack((flux_binned,flux_b))
		
		name="%s_b.lc.txt" % (t[i].split('.')[0])
		z=np.concatenate((time_binned,flux_binned)).reshape(2,len(time_binned)).transpose()
		np.savetxt(name,z,fmt='%.8f    %.3f')
	
		print "[Binning] %d/%d" % ((i+1),len(t))
		
		#pl.plot(time,flux,'r.',time_binned,flux_binned,'g.')
		#pl.show()
		
	return 0

# toggles
bin = 0
bls = 0
analyse_real = 0
analyse_fake = 1
print_results = 0

tdir=["V17.50","V17.75","V18.00","V18.25","V18.50","V18.75","V19.00"]

# bin up lc's to increase accuracy
if bin > 0:
	# run on a few lc's first to test it
	t=cmd.getoutput('ls star_*.lc.txt').split('\n')
	done=BinLcs(t)
	
# edit and run initbls.pro for the star being analyzed
if bls > 0:
	
	for j in range(0,len(tdir)):
		os.chdir(tdir[j])
		os.system("cp /home/jmcc/idl/BLS_Faedi/*.pro .")
		
		t2=cmd.getoutput('ls FakeTransit*.lc.txt').split('\n')
		for i in range(0,len(t2)):
			edit=EditBLSScript(t2[i])
			time.sleep(1)
			os.system('/software/itt/idl71/bin/idl initbls_py')
		
		os.chdir("../")
		
# analyze the real results
if analyse_real > 0:
	t3=cmd.getoutput('ls *_bls.dat').split('\n')
	
	# tokens for final file name
	section = "%s_%s" % (t3[0].split('_')[2],t3[0].split('_')[3])
	
	# arrays for analyses
	sid, Sg, bper, bpow, depth=np.empty(len(t3)),np.empty(len(t3)),np.empty(len(t3)),np.empty(len(t3)),np.empty(len(t3))
	qtran, in1, in2=np.empty(len(t3)),np.empty(len(t3)),np.empty(len(t3))
	
	
	for i in range(0,len(t3)):
		f=open(t3[i]).readlines()
		sid[i]=f[0].split()[0].split('_')[1]
		Sg[i]=f[0].split()[1]
		bper[i]=f[0].split()[2]
		bpow[i]=f[0].split()[3]
		depth[i]=f[0].split()[4]
		qtran[i]=f[1].split()[0]
		in1[i]=f[1].split()[1]
		in2[i]=f[1].split()[2]
	
	name = "BLS_TOTAL_%s.dat" % (section)
	if os.path.exists(name) == True:
		print "Overwritting %s!" % (name)
	
	z=np.concatenate((sid, bper, bpow, depth, qtran)).reshape(5,len(sid)).transpose()
	np.savetxt(name,z,fmt='%05d    %.6f    %.6f    %.6f    %.6f')

	line1=[6]*len(sid)
	
	pl.figure(1)
	pl.plot(sid,Sg,'r-')
	pl.plot(sid,line1,'k-')
	pl.ylim(0,10)
	pl.title('BLS Significance %s' % (section))
	pl.xlabel('Sid')
	pl.ylabel('Significance (Sg)')
	pl.show()
	
	# make a file with those with Sg>6 and their periods for manual phasing. 
	if print_results > 0:
		out_sid,out_Sg,out_bper,out_depth,out_qtran,out_in1,out_in2=[],[],[],[],[],[],[]
		for i in range(0,len(Sg)):
			if Sg[i] > 6.0:
				out_sid.append(sid[i])
				out_Sg.append(Sg[i])
				out_bper.append(bper[i])
				out_depth.append(depth[i])
				out_qtran.append(qtran[i])
				out_in1.append(in1[i])
				out_in2.append(in2[i])
	
		# make numpy copies of array for np.savetxt()
		out_sid=np.copy(out_sid)
		out_Sg=np.copy(out_Sg)
		out_bper=np.copy(out_bper)
		out_depth=np.copy(out_depth)
		out_qtran=np.copy(out_qtran)
		out_in1=np.copy(out_in1)
		out_in2=np.copy(out_in2)
		
		name2 = "BLS_SigDetections_%s.dat" % (section)
		if os.path.exists(name2) == True:
			print "Overwritting %s!" % (name2)
		
		z_fin=np.concatenate((out_sid, out_Sg, out_bper, out_depth, out_qtran, out_in1, out_in2)).reshape(7,len(out_sid)).transpose()
		np.savetxt(name2,z_fin,fmt='%05d    %.6f    %.6f    %.6f    %.6f    %.6f    %.6f')


if analyse_fake > 0:

	for j in range(0,len(tdir)):
		os.chdir(tdir[j])
		t3=cmd.getoutput('ls *_bls.dat').split('\n')
		
		# tokens for final file name
		period=t3[0].split('_')[0][-3:]
		tdur=t3[0].split('_')[1]
		mag = os.getcwd().split('/')[-1]	
			
		# arrays for analyses
		sid, Sg, bper, bpow, depth=np.empty(len(t3)),np.empty(len(t3)),np.empty(len(t3)),np.empty(len(t3)),np.empty(len(t3))
		qtran, in1, in2=np.empty(len(t3)),np.empty(len(t3)),np.empty(len(t3))
		
		
		for i in range(0,len(t3)):
			f=open(t3[i]).readlines()
			sid[i]=f[0].split()[0]
			Sg[i]=f[0].split()[1]
			bper[i]=f[0].split()[2]
			bpow[i]=f[0].split()[3]
			depth[i]=f[0].split()[4]
			qtran[i]=f[1].split()[0]
			in1[i]=f[1].split()[1]
			in2[i]=f[1].split()[2]
		
		x,y,z=0,0,0
		for i in range(0,len(bper)):
			if abs(bper[i]-float(period)) < 0.1:
				if Sg[i] >= 6:
					x=x+1
				if Sg[i] >= 5 and Sg[i] < 6:
					y=y+1
				if Sg[i] < 5:
					z=z+1 
		
		print "[%s] 6 < Sg: %d" % (period,x)
		print "[%s] 5 < Sg < 6: %d" % (period,y)
		print "[%s] Sg < 5: %d" % (period,z)
		print "[%s] Total fraction recovered: %d/%d\n" % (mag,x+y+z,len(bper))
		
		os.chdir("../")

	
	
