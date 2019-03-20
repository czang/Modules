#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator




def get_2_int_lists(gene_file):
	file = open(gene_file,'r')
	List1 = []
	List2 = []
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			List1.append(atoi(sline[0]))
			List2.append(atoi(sline[1]))
	file.close()
	return List1, List2


def add_lists(list1, list2):
	List = []
	assert len(list1) == len(list2)
	for i in range(0,len(list1)):
		t = list1[i] + list2[i]
		List.append(t)
	return List


def main(argv):
	parser = OptionParser()
	parser.add_option("-x", "--inputfile1", action="store", type="string", dest="input_file1", metavar="<file>", help="input file 1")
	parser.add_option("-y", "--inputfile2", action="store", type="string", dest="input_file2", metavar="<file>", help="input file 2")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	list1n, list1x = get_2_int_lists(opt.input_file1)
	list2n, list2x = get_2_int_lists(opt.input_file2)
	assert list1n == list2n
	List = add_lists(list1x, list2x)
	outfile = open(opt.output_file, 'w')
	for i in range(0, len(List)):
		outfile.write(str(list1n[i]) + '\t' + str(List[i]) + '\n')
	outfile.close()


if __name__ == "__main__":
	main(sys.argv)
