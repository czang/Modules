#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

import gene_set_manipulation
import UCSC


def gene_length(file):
	f = open(file, 'r')
	result = {}
	for line in f:
		if not re.match("#",line):
			line = line.strip()
			sline = line.split()
			result[sline[0]] = int(sline[4]) - int(sline[3])
	f.close()
	return result


def get_conversiontable_w_length(infile, length_dic, incol, outcol):
	result = {}
	f = open(infile, 'r')
	for line in f:
		if not re.match("#",line):
			line = line.strip()
			sline = line.split()
			name = sline[incol]
			after = sline[outcol]
			if name in result.keys():
				if length_dic[after] > length_dic[result[name]]:
					result[name] = after
			else:
				result[name] = after
	f.close()
	return result


def convert_id(infile, conversion_table, outfile):
	f = open(infile, 'r')
	o = open(outfile, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			gene_name = sline[0]
			if gene_name in conversion_table.keys():
				sline[0] = conversion_table[gene_name]	
				outline =  '\t '.join(sline) + '\n'
				o.write(outline)
	f.close()
	o.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-k", "--conversion_file", action="store", type="string",
                      dest="conversion_file", help="file with conversion table", metavar="<file>")
	parser.add_option("-g", "--gene_file", action="store", type="string",
                      dest="gene_file", help="input file in bed format", metavar="<file>")
	parser.add_option("-b", "--in_file", action="store", type="string",
                      dest="in_file", help="input file", metavar="<file>")
	parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output file in bed format", metavar="<file>")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
		parser.print_help()
		sys.exit(1)
	
	length_dic = gene_length(opt.gene_file)
	table = get_conversiontable_w_length(opt.conversion_file, length_dic, 1, 0)
	
	f = open(opt.in_file, 'r')
	o = open(opt.outfile, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			old = sline[40]
			if old in table.keys():
				name = table[old]
				o.write(name+'\t'+'\t'.join(sline[:40])+'\n')
			else:
				o.write(old+'\t'+'\t'.join(sline[:40])+'\n')
		else:
			o.write(line)
	f.close()
	o.close()


if __name__ == "__main__":
    main(sys.argv)
