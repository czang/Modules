#!/usr/bin/env python
# Copyright (c) 2012 DFCI & HSPH
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
# chongzhi_zang@dfci.harvard.edu).

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


import gene_set_manipulation
from stats import *


def main(argv):
	parser = OptionParser()
	parser.add_option("-t", "--treatedfile", action="store", type="string", dest="t_file", metavar="<file>", help="treated data file")
	parser.add_option("-b", "--control", action="store", type="int", dest="control", metavar="<int>", help="colum number for name")
	parser.add_option("-f", "--treated", action="store", type="int", dest="treated", metavar="<int>", help="colum number for value")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file")
	
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	
	control_list = gene_set_manipulation.get_float_list(opt.t_file,opt.control)
	treated_list = gene_set_manipulation.get_float_list(opt.t_file,opt.treated)
	treated_names = gene_set_manipulation.get_gene_list(opt.t_file,0)
	control_average = mean(control_list)
	treated_average = mean(treated_list)

	out = open(opt.out_file, 'w')
	for i in range(0,len(control_list)):
		score = sqrt(treated_list[i]/treated_average)-sqrt(control_list[i]/control_average)
		out.write(treated_names[i] + '\t' + str(round(score,5)) + '\n')
	out.close()


if __name__ == "__main__":
	main(sys.argv)
