#!/usr/bin/env python
import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator
import gene_set_manipulation

def categorize(infile, outfile):
	f = open(infile,'r')
	dic = {}
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split('\t')
			assert len(sline) == 2
			domain = sline[0]
			l = sline[1].strip()
			refseq = l.split()
			if domain in dic.keys():
				for item in refseq:
					dic[domain].append(item)
			else:
				dic[domain] = refseq
	f.close()
	o = open(outfile, 'w')
	for domain in dic.keys():
		List = gene_set_manipulation.find_unique_genes(dic[domain])
		for item in List:
			o.write(item + '\t' + domain + '\n')
	o.close()


def main(argv):
	
	parser = OptionParser()
	parser.add_option("-b", "--in_file", action="store", type="string",
                      dest="in_file", help="input file", metavar="<file>")
	parser.add_option("-o", "--output_file", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	categorize(opt.in_file, opt.out_file)


if __name__ == "__main__":
	main(sys.argv)
