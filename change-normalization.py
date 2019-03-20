#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *

infile = open(sys.argv[1]);
norm = atof(sys.argv[2])/atof(sys.argv[3]);


for line in infile:
	line = line.strip();
	sline = line.split();
	outline = sline[0] + "\t" + str(atof(sline[1]) * norm) + "\t" + str(atof(sline[2]) * norm);
	print outline;

infile.close()
