#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *

pc = 1

infile = open(sys.argv[1])

for line in infile:
	if not re.match("#", line):
		line = line.strip()
		sline = line.split()
		L = []
		for i in range(1,len(sline)):
			L.append(atof(sline[i]))
		m = min(L)
		if m <= 0.5:
			for i in range(1,len(sline)):
				s = atof(sline[i]) + pc
				sline[i] = str(s)
			outline = '\t'.join(sline)
		else:
			outline = '\t'.join(sline)
	else:
		outline = line	
	print outline

infile.close()
