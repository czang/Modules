#!/usr/bin/env python
# Copyright (c) 2012 DFCI & HSPH
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
# chongzhi_zang@dfci.harvard.edu).

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

from gene_set_manipulation import *


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--infile", action="store", type="string", dest="infile", metavar="<file>", help="input file")
	parser.add_option("-n", "--groupsize", action="store", type="int", dest="groupsize", metavar="<int>", help="number of elements in one group")
	parser.add_option("-o", "--output_file", action="store", type="string", dest="out_file", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	c = open(opt.infile, 'r')
	total = len(c.readlines())
	c.close()
	groupnumber = total/opt.groupsize
	f = open(opt.infile, 'r')
	g = 1
	i = 0
	while g < groupnumber:
		o = open(opt.out_file+'_'+str(g), 'w')
		for j in range(0, opt.groupsize):
			o.write(f.readline())
			i += 1
		o.close()
		g += 1
	o = open(opt.out_file+'_'+str(g), 'w')
	while i < total:
		o.write(f.readline())
		i +=1
	o.close()
	f.close()


if __name__ == "__main__":
	main(sys.argv)