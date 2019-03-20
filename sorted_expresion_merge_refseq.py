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
	
	
	f = open(opt.in_file, 'r')
	o = open(opt.outfile, 'w')
	head = f.readline()
	if not re.match("#", head):
		head = head.strip()
		firstline = head.split()
		current = firstline[0]
		current_list = []
		n = 1
		for i in range(1,len(firstline)):
			current_list.append(atof(firstline[i]))
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			name = sline[0]
			if name == current:
				for i in range(1,len(sline)):
					current_list[i-1] += atof(sline[i])
				n += 1
			else:
				outlist = []
				for item in current_list:
					outlist.append(str(item/n))
				o.write(current + '\t' + '\t'.join(outlist) + '\n')
				current = name
				current_list = []
				n = 1
				for i in range(1,len(sline)):
					current_list.append(atof(sline[i]))
	outlist = []
	for item in current_list:
		outlist.append(str(item/n))
	o.write(current + '\t' + '\t'.join(outlist) + '\n')
	f.close()
	o.close()


if __name__ == "__main__":
    main(sys.argv)
