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
import scipy
#from Background_island_calculation import *


def factorial(m):
	value = 1.0;
	if m != 0:
		while m != 1:
			value = value*m;
			m = m - 1;
	return value;

def factln(m):
	if m<20:  
		return log(factorial(m));
	else:
		return m*log(m) -m + log(m*(1+4*m*(1+2*m)))/6.0 + log(pi)/2;

def fac(m):
	if m != 0:
		return log(scipy.misc.factorial(m))
	else:
		return -1

def fisher_p(a, b, c, d):
	if a != 0:
		n = a + b + c + d
		return exp(fac(a+b) + fac(c+d) + fac(a+c) + fac(b+d) - fac(a) - fac(b) - fac(c) - fac(d) - fac(n))
	else:
		return -1

def llr(a, k, m, n):
	if a != 0:
		return log(float(a)/float(k)/(float(m)/float(n)), 2)
	else:
		return -100

def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--infile", action="store", type="string", dest="infile", metavar="<file>", help="data file")
	parser.add_option("-m", "--subtotal", action="store", type="int", dest="subtotal", help="subtotal", metavar="<int>")
	parser.add_option("-n", "--total", action="store", type="int", dest="total", help="total", metavar="<int>")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="output file")
	
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	
	m = opt.subtotal
	n = opt.total

	file = open(opt.infile, 'r')
	out = open(opt.out_file, 'w')
	for line in file:
		line = line.strip()
		sline = line.split()
		name = sline[0]
		a = atoi(sline[1])
		k = atoi(sline[2])
		b = k - a
		c = m - a
		d = n - a - b -c
		#p = fisher_p(a, b, c, d)
		l = llr(a, k, m, n)
		#print a, p, l
		out.write(name + '\t' + str(round(l,6)) + '\n')
	file.close()
	out.close()


if __name__ == "__main__":
	main(sys.argv)
