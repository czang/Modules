#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

from gene_set_manipulation import *

def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--genefile1", action="store", type="string", dest="gene_file1", metavar="<file>", help="file with first gene set")
	parser.add_option("-c", "--column-index-for-1", action="store", type="int", dest="c1", metavar="<int>", help="column number for names in the first gene set")
	parser.add_option("-b", "--genefile2", action="store", type="string", dest="gene_file2", metavar="<file>", help="file with second gene set")
	parser.add_option("-u", "--column-index-for-2", action="store", type="int", dest="c2", metavar="<int>", help="column number for names in the second gene set")
	parser.add_option("-o", "--outfilename", action="store", type="string", dest="outfile", metavar="<file>", help="outfile")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)()
	gene_list_1 = get_gene_list(opt.gene_file1, opt.c1);
	gene_list_2 = get_gene_list(opt.gene_file2, opt.c2);
	print len(gene_list_1);
	result = gene_comparison(gene_list_1, gene_list_2);
	
	outfilename1 = opt.outfile + "_shared";
	outfilename2 = opt.outfile + "_disappeared";
	
	ofile1 = open(outfilename1, 'w');
	len_shared = len(result['shared']);
	ofile1.write('# Number of common items between '+ opt.gene_file1 + ' and '+ opt.gene_file2 + " is: " + str(len_shared) + ", out of " + str(len(gene_list_1)) + " items \n");
	for item in result['shared']:	
		ofile1.write(item + '\n');
	ofile1.close();
	
	ofile2 = open(outfilename2, 'w');
	len_disappeared = len(result['not in 2']);
	ofile2.write('# Number of items in '+ opt.gene_file1 + ' but not in '+ opt.gene_file2 + " is: " + str(len_disappeared) + ", out of " + str(len(gene_list_1)) + " items \n");
	for item in result['not in 2']:	
		ofile2.write(item + '\n');
	ofile2.close();
	
	#test = gene_comparison(result['shared'], result['not in 1']);
	#print test['shared'];
	
if __name__ == "__main__":
	main(sys.argv)