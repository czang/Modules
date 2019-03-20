#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


def main(argv):
	parser = OptionParser()
	parser.add_option("-p", "--file1", action="store", type="string", dest="file1", metavar="<file>", help="file 1")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 2:
        	parser.print_help()
        	sys.exit(1)
	
	f = open(opt.file1, 'r')
	h = f.readline()
	h = h.strip()
	sh = h.split()
	length = len(sh) - 1
	f.close()
	m = []
	for i in range(0, length):
		m.append([])
	f = open(opt.file1, 'r')
	for line in f:
		line = line.strip()
		sline = line.split()
		assert len(sline) == length + 1
		for i in range(0, length):
			m[i].append(atof(sline[i+1]))
	f.close()
	for item in m:
		s = sum(item)
		for i in range(0, len(item)):
			item[i] = item[i]/s
	print m


if __name__ == "__main__":
	main(sys.argv)