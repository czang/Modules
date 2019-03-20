#!/usr/bin/env python
# Copyright (c) 2011 DFCI & HSPH
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


import BED
import UCSC
from GenomeData import *
import associate_island_with_genes

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


"input file has to be the output from find-closest-island-for-genes"


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--infile", action="store", type="string", dest="infile", metavar="<file>", help="file from find-closest-island-for-genes")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
		

	f = open(opt.infile, 'r')
	o = open(opt.out_file, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if re.match("None", sline[1]):
				if re.match("None", sline[3]):
					distance = 1e10
				else: 
					distance = int(sline[3])
			elif re.match("None", sline[3]):
				distance = int(sline[1])
			else:
				distance = min(abs(int(sline[1])), abs(int(sline[3])))
			o.write(sline[0]+'\t'+str(distance)+'\n')
	f.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)