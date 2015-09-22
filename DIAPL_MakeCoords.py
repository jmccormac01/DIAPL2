# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_MakeCoords.py - 	a program make .coords files for PSF fitting 
#						
# 						uses tvmark files 
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 28/10/11 - 	code writen
# 28/10/11 -	code tested
#
#	

f=open('M71-tvmark_1_1.coo','r').readlines()

x,y=[],[]
for i in range(0,len(f)):
	x.append(f[i].split()[0])
	y.append(f[i].split()[1])
	
f2=open('M71_1_1.coords','w')

for j in range(0,len(x)):
	line = "%.2f   %.2f\n" % (float(x[j]),float(y[j]))
	f2.write(line)

f2.close()

