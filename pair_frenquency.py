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
from associate_binary_modification_with_expression import *

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def get_pair_set_size(datadic, i, j):
	genelist = datadic.keys()
	resultlist = get_gene_list_with_modification(datadic, genelist, i, 1)
	genelist = resultlist
	resultlist = get_gene_list_with_modification(datadic, genelist, j, 1)
	return len(resultlist)


def item_list(file):
	List = []
	f = open(file, 'r')
	l = f.readline()
	l = l.strip()
	sl = l.split()
	for i in range(1, len(sl)):
		List.append(sl[i])
	f.close()
	return List


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--completedata", action="store", type="string", dest="infile", metavar="<file>", help="input whole data set file name")
	parser.add_option("-n", "--number", action="store", type="int", dest="number", metavar="<int>", help="number of modifications")
	parser.add_option("-o", "--output", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	datadic = read_complete_binary_file(opt.infile)
	items = item_list(opt.infile)
	o = open(opt.outfile, "w")
	for i in range(0, opt.number):
		ni = len(get_gene_list_with_modification(datadic, datadic.keys(), i, 1))
		for j in range(i+1, opt.number):
			nj = len(get_gene_list_with_modification(datadic, datadic.keys(), j, 1))
			np = get_pair_set_size(datadic, i, j)
			print i, j
			o.write(str(items[i]) + '\t' + str(items[j]) + '\t' + str(ni) + '\t' + str(nj) + '\t' + str(np) + '\n')
	o.close()


if __name__ == "__main__":
	main(sys.argv)