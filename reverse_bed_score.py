#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *

infile = open(sys.argv[1]);

for line in infile:
    line = line.strip();
    sline = line.split();
    score = -atof(sline[4]);
    sline[4] = str(score)

    outline = '\t'.join(sline)

    print outline;

infile.close()
