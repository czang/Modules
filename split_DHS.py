#!/usr/bin/env python
# Copyright (c) 2016 DFCI/HSPH
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


def split_regions(infile, outfile, size = 75):
	f = open(infile,'r')
	o = open(outfile, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			start = atoi(sline[1])
			end = atoi(sline[2])
			length = end - start + 1
			flag = (end - start) / size
			if flag <= 1:
				o.write('\t'.join([sline[0], sline[1],sline[2],str(length)])+'\n')
			else:
				l = length / flag
				start_list = []
				start_list.append(start)
				end_list = []
				for i in range(1,flag):
					tag = start + i * l
					end_list.append(tag - 1)
					start_list.append(tag)
				end_list.append(end)
				assert len(start_list) == len(end_list)
				for i in range(0,len(start_list)):
					o.write('\t'.join([sline[0], str(start_list[i]), str(end_list[i]),str(length)])+'\n')

	f.close()
	o.close()


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-b", "--island_bed_file", action="store", type="string",
                      dest="bed_file", help="island bed file", metavar="<file>")
	parser.add_option("-w", "--peakwidth", action="store", type="int", dest="width", help="minimum peak with", metavar="<int>", default=75)
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	split_regions(opt.bed_file, opt.out_file, opt.width)


if __name__ == "__main__":
	main(sys.argv)
