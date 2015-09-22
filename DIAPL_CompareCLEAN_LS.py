

# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_CompareCLEAN_LS.py - 	a python program to compare CLEAN and LS periodograms
#						
# 						
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 07/01/12 - 	code writen
# 07/01/12 -	code tested
#				
#

import numpy as np
import pylab as pl
import commands as cmd

def ComparePlots(output):
	
	# Clean
	freq_c,power_c=np.loadtxt('/Volumes/DATA/NITES/Results/M71/Photometry/Batch1/CLEAN/variables/1_1/star_00544_1_1_wEcS_f2_clean.dat',usecols=[0,1],unpack=True)
	# LS
	freq_ls,power_ls=np.loadtxt('/Volumes/DATA/NITES/Results/M71/Photometry/Batch1/A8.0D12S8C12/1_1/period/star_00544_1_1_scargle.dat', usecols=[0,1],unpack=True)
	
	# periods
	# exlcude the first one because its 0 and makes Inf when 1/0 is ran
	periods_c=1/freq_c[1:]
	periods_ls=1/freq_ls[1:]
	
	# normalise the powers
	# ignore array [0] as before
	power_c=power_c[1:]/(max(power_c))
	power_ls=power_ls[1:]/(max(power_ls))
	
	pl.figure(1)
	pl.subplot(211)
	pl.plot(periods_c,power_c,'r-')
	pl.xlim(0,0.4)
	
	pl.subplot(212)
	pl.plot(periods_ls,power_ls,'k-')
	pl.xlim(0,0.4)
	
	pl.show()
	
	if output > 0:
		name = 'CLEAN_LS_Comp_00544_1_1.lc.txt'
		z=np.concatenate((periods_c[:-1],power_c[:-1],periods_ls,power_ls)).reshape(4,len(periods_ls)).transpose()
		np.savetxt(name,z,fmt='%.5f    %.5f    %.5f    %.5f')

	return 0

compare = 0
output = 0

if compare > 0:
	comp=ComparePlots(output)

# work out difference in periods from CLEAN and LS

#CLEAN
t=cmd.getoutput("ls star*_logfile.dat").split('\n')

# get each CLEAN lc then look for its LS partner

section="%s_%s" % (t[0].split('_')[2],t[0].split('_')[3])
starid=np.empty(len(t))
clean_ps=np.empty(len(t))
ls_ps=np.empty(len(t))

for i in range(0,len(t)):
	starid[i]=t[i].split('_')[1]
	f=open(t[i]).readlines()
	clean_ps[i]=f[6].split()[-1]
	
	f2=open("/Volumes/DATA/NITES/Results/M71/Photometry/Batch1/A8.0D12S8C12/%s/period/star_%05d_%s_logfile.dat" % (section,starid[i],section)).readlines()
	ls_ps[i]=f2[7].split()[-1]
	
diff=clean_ps-ls_ps

for i in range(0,len(diff)):
	print "%05d   %.6f" % (starid[i],diff[i])


