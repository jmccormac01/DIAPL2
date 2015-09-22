

import numpy as np
import pylab as pl
import random 
import commands as cmd
print "Modules loaded!"

# update
# 18/01/2012 - make the transit start time random in the three days before the WF begins


#		3.2d	2.4d	1.6d	1.3Rj		1.6Rj		1.8Rj
#	V	Tdurs1	Tdurs2	Tdurs3	Depths1		Depths2		Depths3		ERR	RMS Plot
# 17.5	2.815	2.531	2.211	0.01075*	0.01628		0.02060		0.007
# 17.75	2.510	2.257	1.972	0.01369*	0.02074		0.02625		0.0085
# 18.0	2.313	2.080	1.817	0.01635*	0.02477		0.03135		0.011
# 18.25	2.139	1.924	1.681	0.01942*	0.02941		0.03722		0.013
# 18.50	1.993	1.792	1.565	0.02280*	0.03454		0.04371		0.017
# 18.75	1.877	1.688	1.474	0.02617*	0.03965		0.05018		0.022
# 19.0	1.783	1.603	1.401	0.02958*	0.04480		0.05671		0.03

######################
# 		toggles		 #
######################
get_wf = 1
make_lcs = 1
plot = 0

Nlcs=50

# set the transit params
#P=3.2
#P=2.4
P=1.6

Depths1=np.array([0.01075,0.01369,0.01635,0.01942,0.02280,0.02617,0.02958])
Depths2=np.array([0.01628,0.02074,0.02477,0.02941,0.03454,0.03965,0.04480])
Depths3=np.array([0.02060,0.02625,0.03135,0.03722,0.04371,0.05018,0.05671])

Tdurs1=np.array([2.815,2.510,2.313,2.139,1.993,1.877,1.783])
Tdurs2=np.array([2.531,2.257,2.080,1.924,1.792,1.688,1.603])
Tdurs3=np.array([2.211,1.972,1.817,1.681,1.565,1.474,1.401])

Err=np.array([0.007,0.0085,0.011,0.013,0.017,0.022,0.03])

Tdur = Tdurs3
Depth = Depths3
ERR = Err

######################

for q in range(0,len(ERR)):
	# N points in transit
	Ntr=int((Tdur[q]/24.0)/0.005)
	if Ntr % 2 != 1:
		Ntr=Ntr + 1
	
	for j in range(0,Nlcs):
		# make time array between start and end times
		# time resolution is 0.005d
		# randomise the starting time to sample the window function properly
		
		tran=random.randrange(5736.0,5739.0)+(random.randrange(0,199)*0.005)
		time=np.arange(tran,5804.005,0.005)
		flux=np.ones(len(time))
		
		# get the time mid points of transit
		diff=(time-time[0]) % P
		n1=np.where(diff<=0.00005)
		
		# work out the number of transit wing points
		if Ntr % 2 == 1:
			wings=(Ntr-1)/2
		if Ntr % 2 == 0:
			wings=Ntr/2
		
		pt1=np.ones(wings-2,float)
		pt2=np.array([0.995,0.98,0.5])
		scale=np.concatenate((pt1,pt2))

		# make the transits
		# be careful at the end if the transit wants to be 
		# made after the array ends 
		i_s=[]
		for i in range(0,len(time)):
			if i==0:
				for k in range(0,wings+1):
					flux[k] = 1 - (scale[k]*Depth[q])
			
			if i != 0 and i in n1[0]:
				for k in range(i,(wings+i)+1):
					if k < len(flux):
						flux[k] = 1 - (scale[k-i]*Depth[q])
				for k in range((i-wings),i):
					if k < len(flux):
						flux[k] = 1 - (scale[-(k-i)]*Depth[q])
		
				i_s.append(i)
		
		#pl.plot(time,flux,'r-')
		#pl.xlim(5735,5755)
		#pl.show()
		
		if get_wf > 0:
			# Get Window Function
			hjd=np.loadtxt('HJDList_2012-01-17.lc.txt',usecols=[0])
			hjd=hjd-2450000
			
			# make and array of 1's which will be set to
			# 0 at times when NITES didn't observe
			time_yn=np.ones(len(time))
			
			diff2=np.empty(len(hjd))
			for i in range(0,len(time)):
				diff2=abs(hjd-time[i])
				if min(diff2) > 0.0001:
					time_yn[i] = 0.0
				print "[%d/%d]" % (i+1,len(time))			
			
			n2=np.where(time_yn==1.0)
			
			if plot > 0:
				pl.plot(time,flux,'r-',time[n2],flux[n2],'g.')
				pl.show()	
			
		
		if make_lcs > 0:
			# add noise
			# noise is got from BLS_FakeTransits.xlsx and level is the average RMS 
			# from the M71 lc's which have RMS < RMSLimit
			flux_new=np.empty(len(flux))
			for i in range(0,len(flux)):
				flux_new[i]=flux[i]+(ERR[q]*random.gauss(0,1.0))
			
			name = 'FakeTransitP%s_Tdur%s_%03d.lc.txt' % (str(P),str(Tdur[q]),j+1)
			f=open(name,'w')
			for i in range(0,len(time[n2])):
				line2 = "%.5f   %.5f\n" % (time[n2][i], flux_new[n2][i])
				f.write(line2)
					
			f.close()
	
				