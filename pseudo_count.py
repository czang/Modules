#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *

infile = open(sys.argv[1])

for line in infile:
	if not re.match("#", line):
		line = line.strip()
		sline = line.split()
		for i in range(1,len(sline)):
			s = round(max(atof(sline[i]),0.02),4)
			sline[i] = str(s)
		outline = '\t'.join(sline)
	else:
		outline = line	
	print outline

infile.close()
