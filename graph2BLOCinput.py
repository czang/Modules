#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import numpy

## get BED module
import BED;
#import UCSC;
import GenomeData;



def chromname(number):
	if number in ch.keys():
		return ch[number]
	else:
		return number


def strandsign(number):
	if number == '1':
		return '+'
	elif number == '-1':
		return '-'
	else:
		return number


def bed2bloc(input_file, chrom, output_file):
	value_list = []
	infile = open(input_file,'r')
	outfile = open(output_file, 'w')
	#gene_list = []
	for line in infile:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			if sline[0] == chrom:
				outfile.write(sline[1] + '\t' + sline[2]  + '\t' + sline[3] + '\n')
				value_list.append(atof(sline[3]))
	infile.close()
	outfile.close()
	return value_list


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfile", action="store", type="string", dest="input_file", metavar="<file>", help="original gene file to be formatted")
	parser.add_option("-b", "--chromosome", action="store", type="string", dest="chrom", metavar="<string>", help="chromosome")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output UCSC format file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	List = bed2bloc(opt.input_file, opt.chrom, opt.output_file)
	Leftover = int(GenomeData.hg18_chrom_lengths[opt.chrom])/1000 + 1 - len(List)
	for i in range(0, Leftover):
		List.append(0.0)
	d = numpy.std(List)
	print "perl ~/BLOC/BLOC_find.pl", opt.output_file, "0", 0.25 * d, ">", opt.chrom+".blocoutput"


if __name__ == "__main__":
	main(sys.argv)
