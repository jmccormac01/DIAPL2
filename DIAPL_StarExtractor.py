# ----------------------------------------------------------------------------------
#								Description
# ----------------------------------------------------------------------------------
#
# DIAPL_StarExtractor.py - 	a program to extract the photometry from split files 
#						
# 						
#

# ----------------------------------------------------------------------------------
# 								Update History
# ----------------------------------------------------------------------------------
# 12/10/11 - 	code writen
# 12/10/11 -	code tested
# 23/11/11 - 	Added flux_err and sky background to star_* files
# 				Changed the names to star_id_section_wES.lc.txt. wES = with ERROR SKY	
#	
#

from numpy import * 
import numpy as np
import commands
import sys
import pyfits

def starExtractor(filePattern, first, nfiles):
	"""This function returns 2 elements: a (nstars x nframes) numpy array storing flux values
	   and a list of nstars elements storing booleans values. In this second list, a "True"
	   means there are null values in the flux somewhere for that star"""
	   
	ret = None  # flux
	ret2 = None # fluxerr
	ret3 = None # sky
	zeroes = None
	
	for j, file in enumerate([(filePattern % x) for x in range(first, first + nfiles)]):
		lines = open(file).readlines()
		if ret is None:
		    ret = np.empty( (len(lines), nfiles) )
		    ret2 = np.empty( (len(lines), nfiles) )
		    ret3 = np.empty( (len(lines), nfiles) )
		    zeroes = [False] * len(lines)

		for i, line in enumerate(lines):
			val1 = line.split()[4]
			val2 = line.split()[5]
			val3 = line.split()[6]
			ret[i, j] = float(val1)
			ret2[i, j] = float(val2)
			ret3[i, j] = float(val3)

			if val1 == "-99999.992" or val1 == 'nan': 
				zeroes[i] = True
		
		print "[%d/%d]" % ((j+1),nfiles)
			
	return ret, ret2, ret3, zeroes


def GetImList():
	f=open('COROT9-imlist-1_1.txt','r')
	s=f.readlines()
	f.close()
	
	imlist=[]
	for i in range(0,len(s)):
		imlist.append(s[i].split('\n')[0])

	return imlist

def GetHJDList(imlist):

	hjdlist=[]
	for i in range(0,len(imlist)):
		
		image=str(imlist[i].split('r')[0]) + ".fits"
		
		path='/data/jmcc/diapl/%s' % (image)	

		hdulist=pyfits.open(path)
		hjdlist.append(hdulist[0].header['HJD-MID'])
	
		print "[%d/%d] HJD got..." % ((i+1),len(imlist))
		
		hdulist.close()
	
	return hjdlist

print "\nProcessing coords files..."
templist=commands.getoutput('ls coords*_1_1.txt').split('\n')
flux, fluxerr, sky, zeroes= starExtractor("coords%05d_1_1.txt", 1, len(templist))

print "\nCoords files processed...\n"

imlist=GetImList()
hjdlist=GetHJDList(imlist)

badstars=[]

for i in range(0,len(flux)):
	
	outflux,outtime,outfluxerr,outsky=[],[],[],[]
	
	for j in range(0,len(templist)):
		if str(flux[i,j]) != "-99999.992" and str(flux[i,j]) != 'nan': 
			outflux.append(flux[i,j])
			outfluxerr.append(fluxerr[i,j])
			outsky.append(sky[i,j])
			outtime.append(hjdlist[j])

	if len(outtime) != len(outflux) != len(outfluxerr) != len(outsky):
		print "Time, Flux, Fluxerr and Sky arrays don't match, exiting!"
		sys.exit()
	
	if len(outtime) >= 100 and len(outflux) >= 100 and len(outfluxerr) >= 100 and len(outsky) >= 100:
	
		name="star_%05d_1_1_wES.lc.txt" % ((i+1))
		file2=open(name,'w')
		
		for k in range(0,len(outtime)):
			line="%s    %s    %s    %s\n" % (outtime[k], outflux[k], outfluxerr[k], outsky[k])
			file2.write(line)
		
		file2.close()
		print "Good star %d/%d" % ((i+1),len(flux))
		
		
	if len(outflux) < 100 and len(outtime) < 100:
		badstars.append((i+1))
		print "Bad star %d/%d [%d]" % ((i+1),len(flux), len(outflux))
		
print "\nTotal stars: %d"  % (len(flux))
print "Bad stars (<100 points): "  + str(len(badstars))		 
	



