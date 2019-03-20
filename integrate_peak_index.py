#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import gene_set_manipulation
import Utility


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--genefile", action="store", type="string", dest="genefile", metavar="<file>", help="gene list file name")
	parser.add_option("-c", "--coeffile", action="store", type="string", dest="infile", metavar="<file>", help="index file prefix")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	List = gene_set_manipulation.get_gene_list(opt.genefile, 0)
	peak_dic = {}
	for gene in List:
		if Utility.fileExists(opt.infile+gene):
			f = open(opt.infile+gene, 'r')
			for line in f:
				if not re.search("coef", line):
					line = line.strip()
					sline = line.split()
					peak = sline[0]
					score = abs(atoi(sline[1]))
					if score < 1000:
						if peak in peak_dic.keys():
							if score < peak_dic[peak]:
								peak_dic[peak] = score
						else:
							peak_dic[peak] = score
			f.close()
	o = open(opt.out_file, 'w')
	for peak in peak_dic.keys():
		o.write(peak + '\t' + str(peak_dic[peak])+'\n')
	o.close()


if __name__ == "__main__":
	main(sys.argv)