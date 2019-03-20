#!/usr/bin/env python
import re, os, sys, shutil
from math import *
from string import *
from optparse import OptionParser
import operator

import BED;
import GenomeData;
import bed_preprocessing;

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
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
		
  	if opt.species not in GenomeData.species_chroms.keys():
		print "The species is not recognized!!";
	else:
		chrom_list = GenomeData.species_chroms[opt.species];
		histogram=[];
		out = open(opt.bed_file +"_multi_copy_tags.bed", 'w');
    		for chrom in chrom_list:
			filename = opt.bed_file + chrom + ".bed";
			if (fileExists(filename)):
				 bed_list = (BED.BED(opt.species, filename, "BED6", -1))[chrom];
				 (plus_bed_list, minus_bed_list) = bed_preprocessing.breakUpStrands(bed_list);
				 print len(plus_bed_list), len(minus_bed_list);
				 plus_bed_list.sort(key=operator.attrgetter('start'));
				 minus_bed_list.sort(key=operator.attrgetter('start'));
				 histogram = bed_preprocessing.combine_histogram(histogram, bed_preprocessing.find_read_copy_distribution(plus_bed_list));
				 histogram = bed_preprocessing.combine_histogram(histogram, bed_preprocessing.find_read_copy_distribution(minus_bed_list));
				 
				 temp_bed_list = bed_preprocessing.find_multi_copy_reads(plus_bed_list, opt.threshold);
				 bed_preprocessing.write_list(temp_bed_list, out);
				 temp_bed_list = bed_preprocessing.find_multi_copy_reads(minus_bed_list, opt.threshold);
				 bed_preprocessing.write_list(temp_bed_list, out);	 
			else:
				print filename + " does not exist";
		bed_preprocessing.write_histogram(histogram, opt.bed_file + "_read_copy_distribution.dat");	
		out.close();
		
if __name__ == "__main__":
	main(sys.argv)
