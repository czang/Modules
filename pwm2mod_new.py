#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


def pwmfile2list(file):
	f = open(file, 'r')
	m = []
	for line in f:
		if not re.search("A",line):
			line = line.strip()
			sline = line.split()
			assert len(sline) == 4
			List = []
			for item in sline:
				List.append(atof(item))
			m.append(List)
	f.close()
	return m


def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--file1", action="store", type="string", dest="file1", metavar="<file>", help="file 1")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	m = pwmfile2list(opt.file1)
	print m


if __name__ == "__main__":
	main(sys.argv)