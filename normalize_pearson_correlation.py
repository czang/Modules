#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import gene_set_manipulation
from GenomeData import *

def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--rawfile", action="store", type="string", dest="raw", metavar="<file>", help="file with raw cross correlation data")
	parser.add_option("-b", "--variancefile1", action="store", type="string", dest="file1", metavar="<file>", help="file with variance data")
	parser.add_option("-c", "--variancefile2", action="store", type="string", dest="file2", metavar="<file>", help="file with the other variance data")
	parser.add_option("-o", "--output", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];

	data_dic = gene_set_manipulation.get_gene_float_dic(opt.raw, 1)
	data_dic_1 = gene_set_manipulation.get_gene_float_dic(opt.file1, 1)
	data_dic_2 = gene_set_manipulation.get_gene_float_dic(opt.file2, 1)
	
	f = open(opt.outfile,'w');
	for chrom in chroms:
		if data_dic.has_key(chrom):
			if data_dic_1.has_key(chrom):
				if data_dic_1[chrom] >= 0.0000001:
					data_dic[chrom] = data_dic[chrom] / sqrt(data_dic_1[chrom])
			if data_dic_2.has_key(chrom):
				if data_dic_2[chrom] >= 0.0000001:
					data_dic[chrom] = data_dic[chrom] / sqrt(data_dic_2[chrom])
			f.write(chrom + '\t' + str(data_dic[chrom]) + '\n')
	f.close();

if __name__ == "__main__":
	main(sys.argv)