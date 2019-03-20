#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser


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
	
	id = 0
	f = open(opt.in_file, 'r')
	o = open(opt.outfile, 'w')
	head = f.readline()
	if not re.match("#", head):
		head = head.strip()
		firstline = head.split()
		current = firstline[id]
		current_value = firstline
		flag = atoi(firstline[4])-atoi(firstline[3])
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			name = sline[id]
			if name == current:
				score = atoi(sline[4])-atoi(sline[3])
				if score > flag:
					current_value = sline
			else:
				o.write('\t'.join(current_value) + '\n')
				current = name
				current_value = sline
				flag = atoi(sline[4])-atoi(sline[3])
	o.write('\t'.join(current_value) + '\n')
	f.close()
	o.close()


if __name__ == "__main__":
    main(sys.argv)
