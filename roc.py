#!/usr/bin/env python
# Copyright (c) 2015 DFCI & HSPH
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


def get_float_list(gene_file, c):
	"""
	c is the 0-based column number 
	Return a list of float numbers
	
	"""
	file = open(gene_file,'r')
	List = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List.append(atof(sline[c]))
	file.close()
	return List


def ROC(List, bin, total):
	x_list = []
	y_list = []
	x_list.append(0.0)
	y_list.append(0.0)
	total_positive = sum(List)
	total_negative = total - total_positive
	call = 0
	true_positive = 0
	false_positive = 0
	for i in range(0,len(List)-1):
		call += bin
		true_positive += List[i]
		false_positive = call - true_positive
		tpr = float(true_positive) / total_positive
		fpr = float(false_positive) / total_negative
		x_list.append(fpr)
		y_list.append(tpr)
	x_list.append(1.0)
	y_list.append(1.0)
	auc = 0.0
	assert len(x_list) == len(y_list)
	for i in range(1,len(x_list)):
		height = (y_list[i-1] + y_list[i]) / 2.0
		width = x_list[i] - x_list[i-1]
		area = height * width
		auc += area
	return x_list, y_list, auc


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input", action="store", type="string", dest="in_file", metavar="<file>", help="input file")
	parser.add_option("-c", "--column", action="store", type="int", dest="column", metavar="<int>", help="column, starting from 1", default=2)
	parser.add_option("-t", "--total", action="store", type="int", dest="total", metavar="<int>", help="total element count")
	parser.add_option("-b", "--binsize", action="store", type="int", dest="bin", metavar="<int>", help="bin size")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="output file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
		parser.print_help()
		sys.exit(1)
	
	inputlist = get_float_list(opt.in_file, opt.column-1)
	xlist, ylist, auc = ROC(inputlist, opt.bin, opt.total)
	print "AUC =", round(auc, 4)
	outfile = open(opt.out_file, 'w')
	for i in range(0,len(xlist)):
		outfile.write(str(round(xlist[i],4)) + '\t' + str(round(ylist[i],4)) + '\n')
	outfile.close()


if __name__ == "__main__":
	main(sys.argv)
