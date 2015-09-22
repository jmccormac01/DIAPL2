# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_Get_rc.py - 	a program to get the cluster radius in arcsec of variables 
#										
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 23/12/11 - 	code writen
# 23/12/11 -	code tested
# 28/05/13 - 	code updated with correct formulae and better way of reading in numbers	
#				formula from here http://www.astronomycafe.net/qadir/q1890.html
#

import numpy as np

# hex to deg
def HexToDeg(x,y):
	if len(x.split(':'))== 3 and len(y.split(':')) == 3:

		ra1,ra2,ra3=x.split(':')
		dec1,dec2,dec3=y.split(':')
		
		ra1=float(ra1)
		ra2=float(ra2)
		ra3=float(ra3)
		
		dec1=float(dec1)
		dec2=float(dec2)
		dec3=float(dec3)
		
		Ra=((ra3/3600.0)+(ra2/60.0)+ra1)*15.0
		Dec=(dec3/3600.0)+(dec2/60.0)+dec1
		
		return Ra, Dec

# toggles
output_rc_file = 1


# read in RA and DEC
f=open('/Users/James/Documents/NITES/M71-SurveyNotes/Positions/AllVariables_positions.cat').readlines()

name,ra,dec=[],[],[]

for i in range(0,len(f)):
	name.append(f[i].split()[0])
	ra.append(f[i].split()[1])
	dec.append(f[i].split()[2])

ra_deg=np.empty(len(ra))
dec_deg=np.empty(len(dec))

# convert to degrees
for i in range(0,len(ra)):
	ra_deg[i],dec_deg[i]=HexToDeg(ra[i],dec[i])
	
	print "%s %s --> %.6f %.6f" % (ra[i],dec[i],ra_deg[i],dec_deg[i])

# cluster centre coords - M71
clus_ra="19:53:46.0"
clus_dec="18:46:40.0"

clus_ra_deg,clus_dec_deg=HexToDeg(clus_ra,clus_dec)

rc=np.empty(len(ra))
rc_1=np.empty(len(ra))
rc_deg=np.empty(len(ra))

for i in range(0,len(ra)):
	rc_1[i]=np.sin(np.radians(clus_dec_deg))*np.sin(np.radians(dec_deg[i])) + np.cos(np.radians(clus_dec_deg))*np.cos(np.radians(dec_deg[i]))*np.cos(np.radians(clus_ra_deg-ra_deg[i]))
	
	rc_deg[i]=np.degrees(np.arccos(rc_1[i]))

	rc[i]=rc_deg[i]*3600.0

	print "%s\t%.10f\t%.6f\t%.6f" % (name[i],rc_1[i],rc_deg[i],rc[i])

if output_rc_file > 0:

	# output the rc values to file
	f2=open('/Users/James/Documents/NITES/M71-SurveyNotes/Positions/AllVariables_rc.txt','w')
	
	f2.write('# ref coords M71 19:53:46.0 +18:46:40.0 - 20130528\n')
	
	for i in range(0,len(rc)):
		line="%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.1f\n" % (name[i],ra[i],dec[i],ra_deg[i],dec_deg[i],rc_deg[i],rc[i])
		
		f2.write(line)
		
	f2.close()
		


