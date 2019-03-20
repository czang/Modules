#!/usr/bin/env python
# Copyright (c) 2014 DFCI/HSPH
# Authors: Chongzhi Zang
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
# chongzhizang@gmail.com).

import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator
import bisect
from gene_set_manipulation import *


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-d", "--densityfile", action="store", type="string",
                      dest="Dfile", metavar="<file>")
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	
	bound = 0.00001
	
	C = get_float_list(opt.Dfile, 2)
	X = get_int_list(opt.Dfile, 0)
	L = []
	i = 0
	while i < len(X):
		l = []
		for j in range(0,50):
			l.append(C[i])
			i += 1
		assert len(l) == 50
		L.append(l)
	assert len(L) == X[-1] - X[0] + 1

	h = []
	for i in range(0,len(L)):
		h.append(L[i][0])
	h_new = [0.0] * (len(L))
	err = 1
	k = 0
	while (err >= bound and k <= 10000):
		E = []
		for i in range(0,len(h)):
			y = 0.0
			for j in range(0,50):
				if i < len(L):
					y += L[i][j]
				if i-50-j >=0:
					y += L[i-50-j][j]
				if i+j+50 < len(h):
					y -= h[i+j+50]
				if i-j-50 >= 0:
					y -= h[i-j-50]
			y = y/50.0
			h_new[i] = y
			E.append(abs(h_new[i]-h[i]))
		err = max(E)
		h = h_new
		k += 1
		print k, err

	
	o = open(opt.out_file, 'w')
	START = X[0]
	for i in range(0, len(h)):
		x = X[0] + i
		y = h[i]
		o.write(str(x) + '\t' + str(y) + '\n')
	print k
	o.close()		


if __name__ == "__main__":
	main(sys.argv)
