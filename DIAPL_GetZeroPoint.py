# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_GetZeroPoint.py - 	a program to get the photometric zero point 
#										
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 21/12/11 - 	code writen
# 21/12/11 -	code tested
#	

import numpy as np

mag=np.array([12.90,13.04,13.47,13.63,13.85,13.90,14.48,14.49,14.54,14.54,14.70])
flux=np.array([283347.9,290177.7,265285.1,288546.3,211352.2,182493.4,96045.77,95728.01,93621.7,96319.07,92834.05])
colour=np.array([0.29,0.58,1.25,1.43,1.42,1.20,1.03,1.05,1.00,1.13,1.21])

zero=np.empty(len(flux))
mag_new=np.empty(len(flux))
alpha=np.empty(len(flux))

for i in range(0,len(flux)):
	zero[i]=mag[i] - (-2.5*np.log10(flux[i]))
	
zp=np.average(zero)

for i in range(0,len(flux)):	
	mag_new[i]=-2.5*np.log10(flux[i])+zp
	alpha[i]=(mag[i]-mag_new[i])/colour[i]
	
alpha_av=np.average(alpha)	

# ZP = 26.974 for batch1 data
print "Zeropoint = %.3f" % (zp)
print "Alpha = %.3f" % (alpha_av)