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

def is_sorted(list):
	"""
	Check if sorted in ascending order.
	input is a list of BED with chrom, start, end and value.
	output: sorted =1 or 0
	"""
	sorted = 1; 
	for index in range(0, len(list)-1):
		if list[index].start > list[index+1].start:
			sorted = 0;
	return sorted; 

def filter_reads(sorted_bed_list, cutoff, outfile):
	"""
	The histogram of tag copy will provide a cutoff for filtering the raw bed, as some of the 
	tags has too many copies.  
	
	return filtered bed objects.
	"""
	if (len(sorted_bed_list) != 0):
		out = open(outfile, 'w')
		#sorted_bed_list.sort(key=operator.attrgetter('start'));
		total_number_tags = len(sorted_bed_list);
		current_value = (sorted_bed_list[0]).start;
		current_count = 1;
		current_tag = sorted_bed_list[0];
		for index in range(1, len(sorted_bed_list)):
			item = sorted_bed_list[index];
			if (item.start != current_value):
				if (current_count <= cutoff):
					write(current_tag, out);
				current_value = item.start;
				current_count = 1; 
				current_tag = item;
			else:
				current_count += 1;
				if (current_count <= cutoff):
					write(current_tag, out);
		if (current_count <= cutoff): #last tag
			write(current_tag, out);		
		out.close();
		
def write(item, out):
	"""
	write one line into outfile. The file openning and closing is handled by outside. 
	"""
	#chrom, start, end, name, score, strand
	outline = item.chrom + "\t" + str(item.start) + "\t" + str(item.end) + "\n";
	out.write(outline);	
	

def main(argv):
	
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="species under consideration", metavar="<str>")
	parser.add_option("-b", "--island_file", action="store", type="string",
                      dest="island_file", help="island file", metavar="<file>")      
	parser.add_option("-o", "--output_file_name", action="store", type="string",
                      dest="out_file", help="output file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
		
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	separate_by_chrom.separateByChrom(chroms, opt.island_file, '.bed1')
	
	for chrom in chroms:
		filename = chrom + ".bed1";
		if (fileExists(filename)):
			bed_list = (BED.BED(opt.species, filename, "BED3", -1))[chrom];
			bed_list.sort(key=operator.attrgetter('start'));
			filter_reads(bed_list, 1, chrom+'.bed2')
	
	separate_by_chrom.combineAllGraphFiles(chroms, '.bed2', opt.out_file)
	separate_by_chrom.cleanup(chroms, '.bed1')
	separate_by_chrom.cleanup(chroms, '.bed2')

if __name__ == "__main__":
	main(sys.argv)
