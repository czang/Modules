#!/usr/bin/env python
# Copyright (c) 2014 DFCI/HSPH
# Authors: Chongzhi Zang
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
# chongzhizang@gmail.com).

import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator




def main(argv):
	
	parser = OptionParser()
	parser.add_option("-i", "--island_file", action="store", type="string",
                      dest="bed_file", help="island bed file", metavar="<file>")
	parser.add_option("-n", "--name", action="store", type="string",
                      dest="name", help="read name", metavar="<string>", default="ID")
	parser.add_option("-c", "--column", action="store", type="int", dest="column", help="column to replace/add the name, 0 based", metavar="<int>", default = 4)
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	c = opt.column
	f = open(opt.bed_file,'r')
	o = open(opt.out_file, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) >= c:
				sline[c] = opt.name
			else:
				sline.append(opt.name)
			o.write('\t'.join(sline) + '\n')
	f.close()
	o.close()


if __name__ == "__main__":
	main(sys.argv)
