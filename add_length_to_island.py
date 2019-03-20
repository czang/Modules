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


def add_tag_density(infile, outfile):
	f = open(infile,'r')
	o = open(outfile, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			start = atoi(sline[1])
			end = atoi(sline[2])
			#count = atof(sline[3])
			length = end - start + 1
			#sline.append(str(length))
			o.write('\t'.join([sline[0], sline[1],sline[2],str(length)])+'\n')
	f.close()
	o.close()


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-b", "--island_bed_file", action="store", type="string",
                      dest="bed_file", help="island bed file", metavar="<file>")
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	add_tag_density(opt.bed_file, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)
