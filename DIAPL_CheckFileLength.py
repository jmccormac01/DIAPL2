# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_CheckFileLength.py - 	a program to check the number of points in a file
#										
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 19/12/11 - 	code writen
# 19/12/11 -	code tested
#
#

import commands as cmd

t=cmd.getoutput('ls star_0*').split('\n')

for i in range(0,len(t)):
	f=open(t[i],'r').readlines()
	
	print "%s:\t%d" % (t[i],len(f))