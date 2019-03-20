#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;
import gene_set_manipulation;


def ratio_calculation(Result_1, Result_2, pseudocount, normalization1, normalization2):
	result = {}
	if len(Result_1.keys()) != len(Result_2.keys()):
		print 'List lengths are not equal!'
	for name in Result_1.keys():
		if Result_2.has_key(name):
			Result_1[name] = Result_1[name] / normalization1 * 10000000
			Result_2[name] = Result_2[name] / normalization2 * 10000000
			if Result_1[name] < 1E-9 or Result_2[name] < 1E-9:
				Result_1[name] += pseudocount
				Result_2[name] += pseudocount
			result[name] = Result_1[name] / Result_2[name]
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--datafile1", action="store", type="string", dest="datafile1", metavar="<file>", help="file with numerator data in the second column")
	parser.add_option("-b", "--datafile2", action="store", type="string", dest="datafile2", metavar="<file>", help="file with denominator data in the second column")
	parser.add_option("-m", "--normalizationfactor1", action="store", type="float", dest="n1", metavar="<float>", help="normalization factor for data 1")
	parser.add_option("-n", "--normalizationfactor2", action="store", type="float", dest="n2", metavar="<float>", help="normalization factor for data 2")
	parser.add_option("-p", "--pseudocount", action="store", type="float", dest="pseudocount", metavar="<float>", help="pseudo count")
	parser.add_option("-l", "--log2", action="store_true", dest="log2", help="take Log2")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
        	parser.print_help()
        	sys.exit(1)
	
	f1 = gene_set_manipulation.get_gene_float_dic(opt.datafile1, 1)
	f2 = gene_set_manipulation.get_gene_float_dic(opt.datafile2, 1)
	
	Result = ratio_calculation(f1, f2, opt.pseudocount, opt.n1, opt.n2)
	
	f = open(opt.outfile,'w')
	for gene in Result.keys():
		r = 0.0
		if opt.log2:
			r = log(Result[gene],2)
		else:
			r = Result[gene]
		f.write(gene + '\t' + str(r) + '\n')
	f.close()


if __name__ == "__main__":
	main(sys.argv)
