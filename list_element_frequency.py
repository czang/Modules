#!/usr/bin/env python
# Copyright (c) 2015 DFCI & HSPH
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


def is_list_sorted(list):
        """
        Check if sorted in ascending order.
        input is a list of values.
        output: sorted =1 or 0
        """
        sorted = 1;
        for index in range(0, len(list)-1):
                if list[index] > list[index+1]:
                        sorted = 0;
        return sorted;


def get_int_list(gene_file, c):
	"""
	c is the 0-based column number 
	Return a list of float numbers
	
	"""
	file = open(gene_file,'r')
	List = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List.append(atoi(sline[c]))
	file.close()
	return List


def frequency(List):
	if not is_list_sorted(List):
		List.sort()
	freq = {}
	i = 0
	while i < len(List):
		flag = List[i]
		c = 1
		i += 1
		while i < len(List):
			if List[i] == flag:
				c += 1
				i += 1
			else:
				freq[flag] = c
				break
		freq[flag] = c
	return freq


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="in_file", metavar="<file>", help="input file")
	parser.add_option("-c", "--column", action="store", type="int", dest="column", metavar="<int>", help="column, starting from 0", default=0)
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="output file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
		parser.print_help()
		sys.exit(1)
	
	inputlist = get_int_list(opt.in_file, opt.column)
	result = frequency(inputlist)
	outfile = open(opt.out_file, 'w')
	for i in result.keys():
		outfile.write(str(i) + '\t' + str(result[i]) + '\n')
	outfile.close()


if __name__ == "__main__":
	main(sys.argv)
