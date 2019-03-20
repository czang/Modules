#!/usr/bin/env python
import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator

import BED;
from GenomeData import *
import bed_preprocessing;
import separate_by_chrom;

Dir = os.getcwd();
grep = "/bin/grep";
cat = "/bin/cat";
plus = re.compile("\+");
minus = re.compile("\-");

def fileExists(f):
	try:
		file = open(f)
	except IOError:
		exists = 0
	else:
		exists = 1
	return exists

def main(argv):
	
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="species under consideration", metavar="<str>")
	parser.add_option("-b", "--raw_bed_file", action="store", type="string",
                      dest="bed_file", help="raw bed file", metavar="<file>")
	parser.add_option("-t", "--threshold", action="store", type="int",
                      dest="threshold", help="threshold for copy number", metavar="<int>")	      
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
		
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	separate_by_chrom.separateByChrom(chroms, opt.bed_file, '.bed1')
	
	for chrom in chroms:
		filename = chrom + ".bed1";
		if (fileExists(filename)):
			bed_list = (BED.BED(opt.species, filename, "BED6", -1))[chrom];
			bed_list.sort(key=operator.attrgetter('start'));
			bed_preprocessing.filter_reads(bed_list, opt.threshold, chrom+'.bed2')
	
	separate_by_chrom.combineAllGraphFiles(chroms, '.bed2', opt.out_file)
	separate_by_chrom.cleanup(chroms, '.bed1')
	separate_by_chrom.cleanup(chroms, '.bed2')

if __name__ == "__main__":
	main(sys.argv)
