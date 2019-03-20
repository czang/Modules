#!/usr/bin/env python
# Copyright (c) 2010 The George Washington University
# Authors: Chongzhi Zang, Weiqun Peng
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
# wpeng@gwu.edu).

import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator




def main(argv):
	
	parser = OptionParser()
	parser.add_option("-i", "--island_file", action="store", type="string",
                      dest="bed_file", help="island bed file", metavar="<file>")
	parser.add_option("-w", "--windowsize", action="store", type="int", dest="window_size", help="bedgraph window size", metavar="<int>", default = 200)
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	w = opt.window_size
	f = open(opt.bed_file,'r')
	o = open(opt.out_file, 'w')
	total = 0
	retained = 0
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			assert len(sline) >= 3
			chrom = sline[0]
			current_start = atoi(sline[1])
			current_end = atoi(sline[2])
			flag = current_start + w
			while flag <= current_end + 1: 
				start = current_start
				end = flag - 1
				o.write(sline[0] + '\t' + str(start) + '\t' + str(end) + '\t1\n')
				current_start = flag
				flag += w
	f.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
