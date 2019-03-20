#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *

infile = open(sys.argv[1]);

for line in infile:
    line = line.strip();
    sline = line.split();
    score = -atof(sline[3]);

    outline = sline[0] + "\t" + sline[1] + "\t" + sline[2] + "\t" + str(score);

    print outline;

infile.close()
