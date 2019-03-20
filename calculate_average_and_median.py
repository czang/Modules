#!/usr/bin/env python
# Copyright (c) 2009 GWU & NHLBI, NIH
# Authors: Chongzhi Zang, Weiqun Peng and Keji Zhao
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
# zang@gwmail.gwu.edu).

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import associate_binary_modification_with_expression


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def get_data_list(infile, colum):
	'''returns a dictionary with geneIDs as keys, expression value as values. colum is the number of colum -1 where the expression data are in the file'''
	file = open(infile, 'r')
	List = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List.append(atof(sline[colum]))
	file.close()
	return List

	
def main(argv):
	parser = OptionParser()
	parser.add_option("-b", "--inputfile", action="store", type="string", dest="bedfile", metavar="<file>", help="input file")
	parser.add_option("-c", "--column", action="store", type="int", dest="column", help="column number, start from 0", metavar="<int>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	List = get_data_list(opt.bedfile, opt.column)
	print 'Number:', len(List), '; Average:', associate_binary_modification_with_expression.average(List), '; Median:',associate_binary_modification_with_expression.middle(List)
	

if __name__ == "__main__":
	main(sys.argv)