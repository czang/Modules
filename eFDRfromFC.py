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
from stats import *
from optparse import OptionParser
import operator
import bisect
import gene_set_manipulation

def eFDRfromFC(L):
	L.sort()
	l = len(L)
	f = 1.85
	while f <= 1.87:
		lgf = log(f,2)
		cut = bisect.bisect_left(L,lgf)
		bg = L[0:cut]
		me = median(bg)
		neg = me * 2 - lgf
		ncut = bisect.bisect_left(L,neg)
		N = l - cut
		fdr = float(ncut)/float(l-cut)
		print f, N, fdr
		f +=0.0005

def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--log2fcfile", action="store", type="string", dest="t_file", metavar="<file>", help="data file")
	
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	
	List = gene_set_manipulation.get_float_list(opt.t_file,1)
	eFDRfromFC(List)

if __name__ == "__main__":
	main(sys.argv)
