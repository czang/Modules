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


def process(infile, outfile):
	f = open(infile,'r')
	o = open(outfile, 'w')
	datalist = []
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split(";")
			sample = {}
			sample['drug'] = sline[0]
			sample['cell'] = sline[1]
			datalist.append(sample)
	o.write('drug\tMCF7\tssMCF7\tHL60\tPC3\tSKMEL5\n')
	current_drug = datalist[0]['drug']
	current_count = [0]*5
	for item in datalist:
		if item['drug'] == current_drug:
			if item['cell'] == "MCF7":
				current_count[0] += 1
			elif item['cell'] == "ssMCF7":
				current_count[1] += 1
			elif item['cell'] == "HL60":
				current_count[2] += 1
			elif item['cell'] == "PC3":
				current_count[3] += 1
			elif item['cell'] == "SKMEL5":
				current_count[4] += 1
		else:
			o.write(current_drug+'\t'+str(current_count[0])+'\t'+str(current_count[1])+'\t'+str(current_count[2])+'\t'+str(current_count[3])+'\t'+str(current_count[4])+'\n')
			current_drug = item['drug']
			current_count = [0]*5
			if item['cell'] == "MCF7":
				current_count[0] += 1
			elif item['cell'] == "ssMCF7":
				current_count[1] += 1
			elif item['cell'] == "HL60":
				current_count[2] += 1
			elif item['cell'] == "PC3":
				current_count[3] += 1
			elif item['cell'] == "SKMEL5":
				current_count[4] += 1
	o.write(current_drug+'\t'+str(current_count[0])+'\t'+str(current_count[1])+'\t'+str(current_count[2])+'\t'+str(current_count[3])+'\t'+str(current_count[4])+'\n')
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
	
	process(opt.bed_file, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)
