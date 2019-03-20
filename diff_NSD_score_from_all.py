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
	parser.add_option("-a", "--controlcenter", action="store", type="string", dest="cc_file", metavar="<file>", help="control center data file")
	parser.add_option("-b", "--treatedcenter", action="store", type="string", dest="tc_file", metavar="<file>", help="treated center data file")
	parser.add_option("-c", "--controlall", action="store", type="string", dest="cf_file", metavar="<file>", help="control all data file")
	parser.add_option("-d", "--treatedall", action="store", type="string", dest="tf_file", metavar="<file>", help="treated all data file")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file")
	
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	
	control_center = gene_set_manipulation.get_float_list(opt.cc_file,1)
	cc_names = gene_set_manipulation.get_gene_list(opt.cc_file,0)
	treated_center = gene_set_manipulation.get_float_list(opt.tc_file,1)
	tc_names = gene_set_manipulation.get_gene_list(opt.tc_file,0)
	assert cc_names == tc_names
	control_flanking = gene_set_manipulation.get_float_list(opt.cf_file,1)
	cf_names = gene_set_manipulation.get_gene_list(opt.cf_file,0)
	treated_flanking = gene_set_manipulation.get_float_list(opt.tf_file,1)
	tf_names = gene_set_manipulation.get_gene_list(opt.tf_file,0)
	assert cf_names == tf_names
	control_average = mean(control_flanking)
	treated_average = mean(treated_flanking)

	out = open(opt.out_file, 'w')
	for i in range(0,len(tf_names)):
		score = sqrt((treated_flanking[i]-treated_center[i])/treated_average)-sqrt((control_flanking[i]-control_center[i])/control_average)-sqrt(treated_center[i]/treated_average)+sqrt(control_center[i]/control_average)
		out.write(cc_names[i] + '\t' + str(round(score,6)) + '\n')
	out.close()


if __name__ == "__main__":
	main(sys.argv)
