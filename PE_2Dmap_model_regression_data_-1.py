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


def float_dic(infile, colum1, colum2):
	file = open(infile, 'r')
	dic = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) >= max(colum1, colum2):
				dic[sline[colum1]] = atof(sline[colum2])
	file.close()
	return dic


def float_list_dic(infile):
	file = open(infile, 'r')
	dic = {}
	List1 = []
	List2 = []
	l = file.readline()
	line = l.strip()
	sline = line.split()
	current = sline[0]
	List1.append(atoi(sline[1]))
	List2.append(atoi(sline[2]))
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			tag = sline[0]
			if tag != current:
				d = {}
				for i in range(0, len(List1)):
					d[List1[i]] = List2[i]
				dic[current] = d
				List1 = []
				List2 = []
				current = tag
				List1.append(atoi(sline[1]))
				List2.append(atoi(sline[2]))
			else:
				assert tag == current
				List1.append(atoi(sline[1]))
				List2.append(atoi(sline[2]))
	d = {}
	for i in range(0, len(List1)):
		d[List1[i]] = List2[i]
	dic[current] = d
	file.close()
	return dic


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-d", "--densityfile", action="store", type="string",
                      dest="Dfile", metavar="<file>")
	parser.add_option("-a", "--h1nsityfile", action="store", type="string",
                      dest="H1file", metavar="<file>")
	parser.add_option("-b", "--h2file", action="store", type="string",
                      dest="H2file", metavar="<file>")
	parser.add_option("-c", "--b1file", action="store", type="string",
                      dest="B1file", metavar="<file>")
	parser.add_option("-f", "--b2file", action="store", type="string",
                      dest="B2file", metavar="<file>")
	parser.add_option("-p", "--Plfile", action="store", type="string",
                      dest="Plfile", metavar="<file>")
	parser.add_option("-e", "--promoterextension", action="store", type="int", dest="promoter_extension", help="upstream and downstream extensions of promoter region", metavar="<int>", default=500)
	parser.add_option("-t", "--TSS", action="store", type="int", dest="TSS", help="TSS", metavar="<int>")
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	
	psuedocount = 0.25
	
	START = opt.TSS - opt.promoter_extension
	END = opt.TSS + opt.promoter_extension - 100
	H1 = float_dic(opt.H1file, 0, 1)
	H2 = float_dic(opt.H2file, 0, 1)
	B1 = float_dic(opt.B1file, 1, 3)
	B2 = float_dic(opt.B2file, 2, 3)
	P = float_dic(opt.Plfile, 0, 1)
	D = float_list_dic(opt.Dfile)
	
	o = open(opt.out_file, 'w')
	
	for i in range(START, END):
		x = str(i)
		h1 = H1[x]
		if h1 >= 2500:
			h1 = log(float(h1)/5000.0, 2)
			b1 = B1[x]
			if x in D.keys():
				d = D[x]
				for j in range(50,100):
					h2 = H2[str(i+j-1)]
					b2 = B2[str(i+j-1)]
					p = P[str(j)]
					if j in d.keys():
						c = float(d[j])
					else:
						c = psuedocount
					#h1 = log(float(h1)/5000.0, 2)
					h2 = log(max(float(h2),1.0)/5000.0, 2)
					p = log(float(p),2)
					c = log(c,2)
					o.write('\t'.join([x, str(j), str(c), str(h1), str(h2), str(b1), str(b2), str(p)]) + '\n')
			else:
				for j in range(50,100):
					h2 = H2[str(i+j-1)]
					h2 = log(max(float(h2),1.0)/5000.0, 2)
					b2 = B2[str(i+j-1)]
					p = P[str(j)]
					p = log(float(p),2)
					c = log(0.25,2)
					o.write('\t'.join([x, str(j), str(c), str(h1), str(h2), str(b1), str(b2), str(p)]) + '\n')
	o.close()		


if __name__ == "__main__":
	main(sys.argv)
