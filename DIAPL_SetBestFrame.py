

# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_SetBestFrame.py	-	A program to set the best frame from template.list
#
#										


# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 18/07/12	- 	code writen
# 18/07/12	-	code tested
#

s=open('template.list').readlines()
t=open('diapl_setup.par').readlines()

refim=s[0].split()[0]
t[41]='REFIM="%s"\n' % (refim)

f=open('diapl_setup.par','w')

for i in range(0,len(t)):
	f.write(t[i])

f.close()
