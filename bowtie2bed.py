#!/usr/bin/env python
# Copyright (c) 2010 DFCI/HSPH
# Authors: Chongzhi Zang and X. Shirley Liu
#
# This software is distributable under the terms of the GNU
# General Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or otherwise
# using this module constitutes acceptance of the terms of this License.
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
# 

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;



"""
This module is used to convert a Bowtie output format file in to a BED file.

"""

plus = re.compile('\+');
minus = re.compile('\-');

def chromname(number):
	if number in ch.keys():
		return ch[number]
	else:
		return number


def strandsign(number):
	if number == '1':
		return '+'
	elif number == '-1':
		return '-'
	else:
		return number

def bowtie2BED(input_file, output_file):
	testinfile = open(input_file,'r')
	# Figure out which column to start
	test = testinfile.readline()
	test = test.strip()
	test = test.split()
	i = 0
	while i < len(test):
		if plus.match(test[i]) or minus.match(test[i]):
			strand = i
			break
		else:
			i += 1
	chrom = strand + 1
	start = strand + 2
	seq = strand + 3
	rest = strand + 5
	tagsize = len(test[seq])
	testinfile.close()
	
	infile = open(input_file,'r')
	outfile = open(output_file, 'w')
	for line in infile:
		line = line.strip()
		sline = line.split()
		if plus.match(sline[strand]):
			left = sline[start]
			right = str(int(sline[start]) + tagsize)
		elif minus.match(sline[strand]):
			left = str(int(sline[start]) - tagsize)
			right = sline[start]
		outfile.write(sline[chrom] + '\t' + left  + '\t' + right + '\t' + sline[seq] + '\t' + sline[rest] + '\t' + sline[strand] + '\n')
	infile.close()
	outfile.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="bowtie output format file")
#	parser.add_option("-f", "--tagsize", action="store", type="int", dest="tagsize", metavar="<int>", help="sequence read length", default=25)
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output BED file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	bowtie2BED(opt.input_file, opt.output_file)


if __name__ == "__main__":
	main(sys.argv)
