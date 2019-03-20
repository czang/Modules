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
from scipy import *
from string import *
from optparse import OptionParser
import operator

#import associate_binary_modification_with_expression


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def binomial_probability_calculator(p, n, k):
	return (round(misc.comb(n, k)) * pow(p, k) * pow(1 - p, n - k))


def main(argv):
	parser = OptionParser()
	parser.add_option("-m", "--total", action="store", type="int", dest="m", help="pool size, natural number", metavar="<int>")
	parser.add_option("-p", "--total_positive", action="store", type="int", dest="p", help="total positive number to decide p, must be less than m", metavar="<int>")
	parser.add_option("-n", "--n", action="store", type="int", dest="n", help="n", metavar="<int>")
	parser.add_option("-k", "--k", action="store", type="int", dest="k", help="k", metavar="<int>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	print 'Binomial probability is', binomial_probability_calculator(float(opt.p)/opt.m, opt.n, opt.k)
	

if __name__ == "__main__":
	main(sys.argv)