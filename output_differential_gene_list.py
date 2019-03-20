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


## from a gene expression data matrix, output a list of genes with at least one sample abs >= threshold. 


def process(infile, threshold, outfile):
	f = open(infile,'r')
	o = open(outfile, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List = []
			for i in range(1,len(sline)):
				List.append(abs(atof(sline[i])))
			if max(List) >= threshold:
				o.write(sline[0]+'\n')
	f.close()
	o.close()


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-a", "--all_input_file", action="store", type="string",
                      dest="in_file", help="island bed file", metavar="<file>")
	parser.add_option("-t", "--threshold", action="store", type="float",
                      dest="threshold", help="cutoff", metavar="<float>")
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	process(opt.in_file, opt.threshold, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)
