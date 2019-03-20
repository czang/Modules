#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import scipy.stats

## get BED module
import BED;
import UCSC;
import gene_set_manipulation;
import islands_statistics_pr;
import associate_binary_modification_with_expression;

def pdf2cdf(infile, outfile):
	x = []
	y = []
	f=open(infile, 'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if len(sline) >= 2:
				x.append(int(sline[0]))
				y.append(float(sline[1]))
	f.close()
	total = sum(y)
	cum = []
	for item in y:
		total = total - item
		cum.append(total)
	assert len(x) == len(cum)
	o = open(outfile, 'w')
	o.write('# CDF of island length distribution\n')
	o.write('0\t' + str(sum(y)) + '\n')
	for i in range(0, len(x)):
		o.write(str(x[i]) + '\t' + str(cum[i]) + '\n')
	o.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original gene file to be formatted")
	parser.add_option("-t", "--type", action="store", type="int", dest="option", help="1 for PDF, 2 for CDF", metavar="<int>")     
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output UCSC format file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6 or opt.option not in [1,2]:
        	parser.print_help()
        	sys.exit(1)
	if opt.option == 1: 
		islands_statistics_pr.find_islands_length_histogram({}, opt.input_file, opt.output_file)
	elif opt.option == 2: 
		islands_statistics_pr.find_islands_length_histogram({}, opt.input_file, "pdftemp")
		pdf2cdf("pdftemp", opt.output_file)
		try:
			if os.remove('%s' %("pdftemp")): raise
		except: sys.stderr.write("clean up failed\n");
	
	length_list = []
	infile = open(opt.input_file, 'r');
	for line in infile:
		if not re.match("track", line):
			line = line.strip();
			sline = line.split();
			assert (len(sline) >= 3);
			length = atoi(sline[2]) - atoi(sline[1]) + 1;
			assert (length >= 0);
			length_list.append(length)
	infile.close();
	#print "Total number of islands:", len(length_list)
	#print "Total length of islands:", sum(length_list)
	#print "Average length of islands:", scipy.stats.mean(length_list)
	#print "Median length of islands:", scipy.stats.median(length_list)
	print opt.input_file, "Ntotal", len(length_list), "Ltotal", sum(length_list), "Lavg", scipy.stats.mean(length_list), "Lmed", scipy.stats.median(length_list)



if __name__ == "__main__":
	main(sys.argv)
