#!/usr/bin/env python
# Copyright (c) 2012 DFCI & HSPH
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

from gene_set_manipulation import *


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--wiggle", action="store", type="string", dest="in_file", metavar="<file>", help="input file")
	parser.add_option("-c", "--target", action="store", type="string", dest="targetlist", metavar="<file>", help="target list file name")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="output file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	inputlist = get_gene_list(opt.in_file, 0)
	targets = get_gene_list(opt.targetlist, 0)
	outfile = open(opt.out_file, 'w')
	for item in inputlist:
		if item in targets:
			outfile.write(item + '\t1\n')
		else:
			outfile.write(item + '\t0\n')
	outfile.close()


if __name__ == "__main__":
	main(sys.argv)
