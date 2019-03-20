#!/usr/bin/env python
# Copyright (c) 2010 DFCI/HSPH
# Author: Chongzhi Zang
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


def process(file1, file2, outfile):
	f = open(file1,'r')
	g = open(file2, 'r')
	o = open(outfile, 'w')
	dic1 = {}
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			dic1[sline[0]] = atoi(sline[1])
	for line in g:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			id = sline[0]
			score = atoi(sline[1])
			o.write(id+'\t'+string(dic1[id]+score)+'\n')
	f.close()
	g.close()
	o.close()


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-a", "--file1", action="store", type="string",
                      dest="file1", help="file1", metavar="<file>")
	parser.add_option("-b", "--file2", action="store", type="string",
                      dest="file2", help="file2", metavar="<file>")
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	process(opt.file1, opt.file2, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)
