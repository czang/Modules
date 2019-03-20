#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser


def main(argv):
	parser = OptionParser()
	parser.add_option("-g", "--gene_file", action="store", type="string",
                      dest="in_file", help="input file", metavar="<file>")
	parser.add_option("-l", "--refseqID", action="store", type="int", dest="id", metavar="<int>", help="colum number for ID, start from 0, must be sorted")
	parser.add_option("-f", "--value", action="store", type="int", dest="value", metavar="<int>", help="colum number for value, start from 0")
	parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output file", metavar="<file>")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
		parser.print_help()
		sys.exit(1)
	
	
	f = open(opt.in_file, 'r')
	o = open(opt.outfile, 'w')
	head = f.readline()
	if not re.match("#", head):
		head = head.strip()
		firstline = head.split()
		current = firstline[opt.id]
		current_value = firstline
		flag = atof(firstline[opt.value])
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			name = sline[opt.id]
			if name == current:
				score = atof(sline[opt.value])
				if score > flag:
					current_value = sline
			else:
				o.write('\t'.join(current_value) + '\n')
				current = name
				current_value = sline
				flag = atof(sline[opt.value])
	o.write('\t'.join(current_value) + '\n')
	f.close()
	o.close()


if __name__ == "__main__":
    main(sys.argv)
