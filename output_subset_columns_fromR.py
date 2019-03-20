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


## from a input data list of id, concentration, drug name from one cell type. output 1019 drug samples: id, drug. 


def get_lists(infile):
	f = open(infile, 'r')
	List1 = []
	List2 = []
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List1.append(sline[0])
			List2.append(sline[1])
	f.close()
	return List1, List2


def process(infile, id_list, label_list, outfile):
	f = open(infile,'r')
	o = open(outfile, 'w')
	topline = f.readline()
	topline = topline.strip()
	sline = topline.split()
	L = len(sline)
	index_list = []
	index_list.append(0)
	for item in id_list:
		id = 'X'+item
		if id in sline:
			index_list.append(sline.index(id))
		else:
			print "error!"
	o.write(' \t'+'\t'.join(label_list)+'\n')
	for line in f:
		line = line.strip()
		sline = line.split()
		assert len(sline) == L
		for index in index_list:
			o.write(sline[index]+'\t')
		o.write('\n')
	f.close()
	o.close()


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-a", "--all_input_file", action="store", type="string",
                      dest="in_file", help="island bed file", metavar="<file>")
	parser.add_option("-b", "--druglist_file", action="store", type="string",
                      dest="drug_file", help="island bed file", metavar="<file>")
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	(id_list, label_list) = get_lists(opt.drug_file)
	#print len(druglist)
	process(opt.in_file, id_list, label_list, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)
