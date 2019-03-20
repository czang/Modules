#!/usr/bin/env python
# Copyright (c) 2010 DFCI/HSPH
# Authors: Chongzhi Zang and Xiaole Shirley Liu
#
# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
#
# Disclaimer
# 
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# czang@jimmy.harvard.edu).

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect


Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--file", action="store", type="string", dest="in_file", metavar="<file>", help="input file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	f = open(opt.in_file, 'r')
	xnow = 0
	ynow = 0.0
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			x = atoi(sline[0])
			y = atof(sline[1])
			if y > ynow: 
				xnow = x
				ynow = y
	f.close()
	print opt.in_file, xnow
			

if __name__ == "__main__":
	main(sys.argv)