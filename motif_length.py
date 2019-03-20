#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *

cutoff = 15

infile = open(sys.argv[1])

line = infile.readline()
line = line.strip()
sline = line.split()
size = atoi(sline[2])-atoi(sline[1])
if size >= cutoff:
	print sys.argv[1]
infile.close()
