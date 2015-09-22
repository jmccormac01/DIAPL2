# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_AnalyseSignificance.py - 	a python program to work on significance_section.txt
#									files
#
#		This is primarily to look for SX PHe type objects. 
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 06/12/11 - 	code writen
# 06/12/11 -	code tested
#				
#

import commands as cmd
import numpy as np

t=cmd.getoutput('ls significance*.txt').split('\n')

section="%s_%s" % (t[0].split("_")[1],t[0].split("_")[2])

starid,peak_period,peak_power,Sg,Sg_fac=np.loadtxt(t[0],usecols=[0,1,2,3,4],unpack=True)

gt16=[]
gt10=[]
gt4=[]

for i in range(0,len(Sg_fac)):
	if Sg_fac[i] >= 16:
		gt16.append(i)
	
	if Sg_fac[i] >= 10 and Sg_fac[i] < 16:
		gt10.append(i)

	if Sg_fac[i] >= 4 and Sg_fac[i] < 10:
		gt4.append(i)

print "Stars Sg_fac >= 16:"
for i in range(0,len(gt16)):
	print "Star: %d Peak Period: %.6f Sg_fac: %.3f" % (starid[gt16[i]],peak_period[gt16[i]],Sg_fac[gt16[i]])
	
print "Stars 10 <= Sg_fac < 16 :"
for i in range(0,len(gt10)):
	print "Star: %d Peak Period: %.6f Sg_fac: %.3f" % (starid[gt10[i]],peak_period[gt10[i]],Sg_fac[gt10[i]])	
	
print "Stars 4 <= Sg_fac < 10:"
for i in range(0,len(gt4)):
	print "Star: %d Peak Period: %.6f Sg_fac: %.3f" % (starid[gt4[i]],peak_period[gt4[i]],Sg_fac[gt4[i]])	
	
	
	
	
	
	
	
	
	
	

		
		