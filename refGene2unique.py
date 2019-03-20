#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

import gene_set_manipulation
import UCSC


def main(argv):
	parser = OptionParser()
	parser.add_option("-g", "--gene_file", action="store", type="string",
                      dest="in_file", help="input file in bed format", metavar="<file>")
	parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output file in bed format", metavar="<file>")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
		parser.print_help()
		sys.exit(1)
	
	
	f = open(opt.in_file, 'r')
	length_dic = {}
	line_dic = {}
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			name = sline[1]
			length = int(sline[5]) - int(sline[4])
			if not name in length_dic.keys():
				length_dic[name] = length
				line_dic[name] = '\t'.join(sline[1:])
			elif length > length_dic[name]:
				length_dic[name] = length
				line_dic[name] = '\t'.join(sline[1:])
		else:
			line = line.strip()
			sline = line.split()
			headline = "#" + '\t'.join(sline[1:])
	f.close()
	o = open(opt.outfile, 'w')
	o.write(headline + '\n')
	for name in line_dic.keys():
		o.write(line_dic[name] + '\n')
	o.close()


if __name__ == "__main__":
    main(sys.argv)
