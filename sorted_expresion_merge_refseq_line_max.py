#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser


def list_max(List):
	# List is a list of string formatted numbers, return the maximum number
	L = []
	for item in List:
		L.append(atof(item))
	return max(L)


def main(argv):
	parser = OptionParser()
	parser.add_option("-g", "--gene_file", action="store", type="string",
                      dest="in_file", help="input file", metavar="<file>")
	parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output file", metavar="<file>")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
		parser.print_help()
		sys.exit(1)
	
	
	f = open(opt.in_file, 'r')
	o = open(opt.outfile, 'w')
	head = f.readline()
	if not re.match("#", head):
		head = head.strip()
		firstline = head.split()
		current = firstline[0]
		current_list = firstline[1:len(firstline)]
		flag = list_max(current_list)
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			name = sline[0]
			if name == current:
				score = list_max(sline[1:len(sline)])
				if score > flag:
					current_list = sline[1:len(sline)]
					flag = score
			else:
				o.write(current + '\t' + '\t'.join(current_list) + '\n')
				current = name
				current_list = sline[1:len(sline)]
				flag = list_max(current_list)
	f.close()
	o.close()


if __name__ == "__main__":
    main(sys.argv)
