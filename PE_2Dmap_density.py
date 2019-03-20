#!/usr/bin/env python
# Copyright (c) 2014 
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

import GenomeData
import SeparateByChrom
import Utility


def unique_density_sorted(infile, outfile):
	'''infile can only contain reads from one chromosome and only one kind of strands (+/-). file must be pre-sorted by column2, then by column3.'''
	f = open(infile,'r')
	o = open(outfile, 'w')
	total = 0
	retained = 0
	l = f.readline()
	l = l.strip()
	sline = l.split()
	if len(sline) >= 3:
		current_x = atoi(sline[1])
		current_y = atoi(sline[3])
		current_count = 1
		total = 1
		retained = 0
		for line in f:
			if not re.match("#", line):
				total += 1
				line = line.strip()
				sline = line.split()
				x = atoi(sline[1])
				y = atoi(sline[3])
				if x != current_x:
					o.write('\t'.join([str(current_x),str(current_y),str(current_count)]) + '\n')
					retained += 1
					current_line = sline
					current_x = x
					current_y = y
					current_count = 1
				elif y != current_y:
					o.write('\t'.join([str(current_x),str(current_y),str(current_count)]) + '\n')
					retained += 1
					current_line = sline
					current_x = x
					current_y = y
					current_count = 1
				else:
					current_count += 1
					assert current_x == x
					assert current_y == y
		o.write('\t'.join([str(current_x),str(current_y),str(current_count)]) + '\n')
		retained += 1
	f.close()
	o.close()
	return (total, retained)


def process(infile, outfile):
	'''infile can only contain reads for one gene'''
	try:
		if os.system('sort -g -k 2,4 %s > temp' % (infile)):
			raise
	except: sys.stderr.write(str(infile) + " reads do not exist\n");
	(p_total, p_retained) = unique_density_sorted('temp', outfile)
	#print chrom, "\tAll fragments:",p_total, "\tRetained fragments:", p_retained 
	os.system('rm temp')


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-b", "--raw_bed_file", action="store", type="string",
                      dest="bed_file", help="raw bed file", metavar="<file>")
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	process(opt.bed_file, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)
