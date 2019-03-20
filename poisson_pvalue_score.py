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


import scipy.stats
import gene_set_manipulation


def poisson_pvalue_score(observation, average):
	if observation < 0.2:
		observation = 0.2
	if average < 0.2:
		average = 0.2
	if observation > average:
		pvalue = scipy.stats.poisson.sf(observation, average)[()]
		return -log10(pvalue)
	elif observation < average:
		pvalue = scipy.stats.poisson.sf(average, observation)[()]
		return log10(pvalue)
	else:
		return 0.0


def main(argv):
	parser = OptionParser()
	parser.add_option("-c", "--controlfile", action="store", type="string", dest="c_file", metavar="<file>", help="control data file")
	parser.add_option("-t", "--treatedfile", action="store", type="string", dest="t_file", metavar="<file>", help="treated data file")
	parser.add_option("-a", "--controlsize", action="store", type="float", dest="control_size", metavar="<float>", help="control sequencing depth")
	parser.add_option("-b", "--treatedsize", action="store", type="float", dest="treated_size", metavar="<float>", help="treated sequencing depth")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file")
	
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	
	control_list = gene_set_manipulation.get_float_list(opt.c_file,1)
	control_names = gene_set_manipulation.get_gene_list(opt.c_file,0)
	treated_list = gene_set_manipulation.get_float_list(opt.t_file,1)
	treated_names = gene_set_manipulation.get_gene_list(opt.t_file,0)
	scaling_factor = opt.treated_size*1.0/opt.control_size
	assert control_names == treated_names

	out = open(opt.out_file, 'w')
	for i in range(0,len(control_list)):
		observation = treated_list[i]
		average = control_list[i] * scaling_factor
		score = poisson_pvalue_score(observation, average)
		out.write(control_names[i] + '\t' + str(round(score,5)) + '\n')
	out.close()


if __name__ == "__main__":
	main(sys.argv)
