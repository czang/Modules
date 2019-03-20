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
	parser.add_option("-f", "--treated", action="store", type="int", dest="treated", metavar="<int>", help="colum number for value")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file")
	
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	
	control_list_1 = gene_set_manipulation.get_float_list(opt.t_file,1)
	control_list_2 = gene_set_manipulation.get_float_list(opt.t_file,2)
	control_list_3 = gene_set_manipulation.get_float_list(opt.t_file,3)
	treated_list = gene_set_manipulation.get_float_list(opt.t_file,opt.treated)
	treated_names = gene_set_manipulation.get_gene_list(opt.t_file,0)
	control_average_1 = mean(control_list_1)
	control_average_2 = mean(control_list_2)
	control_average_3 = mean(control_list_3)
	treated_average = mean(treated_list)

	out = open(opt.out_file, 'w')
	for i in range(0,len(control_list_1)):
		score1 = sqrt(treated_list[i]/treated_average)-sqrt(control_list_1[i]/control_average_1)
		score2 = sqrt(treated_list[i]/treated_average)-sqrt(control_list_2[i]/control_average_2)
		score3 = sqrt(treated_list[i]/treated_average)-sqrt(control_list_3[i]/control_average_3)
		score = (score1 + score2 + score3) / 3
		out.write(treated_names[i] + '\t' + str(round(score,5)) + '\n')
	out.close()


if __name__ == "__main__":
	main(sys.argv)
