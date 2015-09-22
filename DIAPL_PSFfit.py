# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# PSF.py - 	a program to extract the magnitudes from the template PSF fitted image
#						
# 						
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 17/10/11 - 	code writen
# 17/10/11 -	code tested
#
#	
	

from pylab import *

f=open('tplr1_1.als.1','r')
s=f.readlines()

num,x,y,mag,magerr=[],[],[],[],[]

s=s[44:]

for i in range(0,len(s)/2):
	num.append(int(s[::2][i].split()[0]))
	x.append(s[::2][i].split()[1])
	y.append(s[::2][i].split()[2])
	if s[::2][i].split()[3] == 'INDEF':
		mag.append(1.0)
	if s[::2][i].split()[3] != 'INDEF':	
		mag.append(float(s[::2][i].split()[3]))
		magerr.append(float(s[::2][i].split()[4]))
		
sorted=zip(num,x,y,mag,magerr)
sorted.sort()

num2,x2,y2,mag2,magerr2=[],[],[],[],[]

for i in range(0,len(sorted)):
	num2.append(sorted[i][0])
	x2.append(sorted[i][1])
	y2.append(sorted[i][2])
	mag2.append(sorted[i][3])
	magerr2.append(sorted[i][4])
	
flux=[]

for i in range(0,len(sorted)):
	x=(mag2[i]-25.0)/(-2.5)
	flux.append(pow(10.0,x))
	