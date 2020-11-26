#!/usr/local/bin/python2.6
import sys
import math
import os
import operator
from numpy import NaN, Inf, arange, isscalar, asarray

#### add an extra chr and at least 2-3 data points to the end of the wig file (these should be deleted afterwards from the file it generates)
#example command:
#python peakdetect.py lactis_0.985.wig 0.1 > peaks_lactis_0.985_0.1.txt (the wig file is produced if you put the smoothing on in our plotting script; the number states the minimal peak height that should be considered)



def peakdet(v, delta, x = None):
    """
Converted from MATLAB script at http://billauer.co.il/peakdet.html
Currently returns two lists of tuples, but maybe arrays would be better
function [maxtab, mintab]=peakdet(v, delta, x)
%PEAKDET Detect peaks in a vector
% [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
% maxima and minima ("peaks") in the vector V.
% MAXTAB and MINTAB consists of two columns. Column 1
% contains indices in V, and column 2 the found values.
%
% With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
% in MAXTAB and MINTAB are replaced with the corresponding
% X-values.
%
% A point is considered a maximum peak if it has the maximal
% value, and was preceded (to the left) by a value lower by
% DELTA.
% Eli Billauer, 3.4.05 (Explicitly not copyrighted).
% This function is released to the public domain; Any use is allowed.
"""
    maxtab = []
    mintab = []
       
    if x is None:
        x = arange(len(v))
    
    v = asarray(v)
    
    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')
    
    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')
    
    if delta <= 0:
        sys.exit('Input argument delta must be positive')
    
    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN
    
    lookformax = True
    
    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]
        
        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return maxtab, mintab

filename = sys.argv[1]
delta = sys.argv[2]
datafile = open(filename, 'r')

data=datafile.readlines()
datafile.close()
x=[]
v=[]


header = ''
sol = ''
for line in data:
	if line.find('variableStep') != -1:
		if header != '':
			print header
			sol = peakdet(v,float(delta),x)
    		
			for item in sol[0]:
				print str(item[0]) + "\t" + str(item[1])
    			
		header = line.strip()
		x = []
		v = []
	else:
		val=line.split()
		x.append(float(val[0]))
		v.append(float(val[1]))

print header
for item in sol[0]:
	print str(item[0]) + "\t" + str(item[1])
