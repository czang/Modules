#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *

infile = open(sys.argv[1]);

for line in infile:
	line = line.strip();
	sline = line.split();
	#score = -atoi(sline[3]);
	if sline[2] == 'F':
		outline = sline[0] + "\t" + sline[1] + "\t" + str(atoi(sline[1])+25) + "\t0\t0\t+";
	elif sline[2] == 'R':
		outline = sline[0] + "\t" + str(atoi(sline[1])-25) + "\t" + sline[1] + "\t0\t0\t-";
	print outline;

infile.close()
