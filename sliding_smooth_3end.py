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
import bisect


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-i", "--bedfile", action="store", type="string",
                      dest="bed_file", help="island bed file", metavar="<file>")
	parser.add_option("-t", "--tss", action="store", type="int", dest="TSS", help="TSS coordinate", metavar="<int>")
	parser.add_option("-e", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream extensions of promoter region", metavar="<int>", default=800)
	parser.add_option("-w", "--windowsize", action="store", type="int", dest="window", help="windowsize", metavar="<int>", default=50)
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	ext = opt.window/2
	START = opt.TSS - opt.promoter_extension
	END = opt.TSS + opt.promoter_extension
	f = open(opt.bed_file,'r')
	o = open(opt.out_file, 'w')
	read_list = []
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			read_list.append(atoi(sline[2]))
	f.close()
	read_list.sort()
	for i in range(START, END):
		i_start = i - ext
		i_end = i + ext
		s = bisect.bisect_left(read_list, i_start)
		e = bisect.bisect_right(read_list, i_end)
		count = e - s
		o.write(str(i) + '\t' + str(count) + '\n')
	o.close()


if __name__ == "__main__":
	main(sys.argv)
