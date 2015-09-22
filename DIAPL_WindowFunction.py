# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_WindowFunction.py - 	a python program to plot window function
#
#		
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 13/01/12 - 	code writen
# 13/01/12 -	code tested
#				
#

import pylab as pl
import numpy as np

output = 1

freq,power=np.loadtxt('/Volumes/DATA/NITES/Results/M71/Photometry/Batch1/A6.0D12S8C12/1_1/wES/filtered2/windowfunction.dat',usecols=[0,1],unpack=True)

# trim the arrys to P<10d
freq_low=np.copy(freq[40:])
power_low=np.copy(power[40:])
period_low=1/freq_low 

# normalise the powers to 1
power_low=power_low/(max(power_low))

pl.plot(period_low,power_low,'r-')
pl.plot(period_low,power_low,'r-')
pl.show()

if output > 0:
	name = "windowfunction_2012-01-13.lc.txt"
	z=np.concatenate((period_low,power_low)).reshape(2,len(period_low)).transpose()
	np.savetxt(name,z,fmt="%.5f    %.5f")
	