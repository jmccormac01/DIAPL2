

# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_SetImageList.py	-	A program to set the image list and
#							FIELD in diapl_setup.par for current data
#
#										


# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 22/07/12	- 	code writen
# 22/07/12	-	code tested
#

import commands as cmd

t=open('diapl_setup.par').readlines()

FIELD=cmd.getoutput('ls *-*-*.fits').split('\n')[0].split('-')[0]
IMAGES=cmd.getoutput('ls *-*-Images.txt').split('\n')[0]

t[34]='IMAGES="%s"\n' % (IMAGES)
t[40]='FIELD="%s"\n' % (FIELD)

f=open('diapl_setup.par','w')

for i in range(0,len(t)):
	f.write(t[i])

f.close()
