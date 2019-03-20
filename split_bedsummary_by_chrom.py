#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

#from species_chroms import *
import GenomeData;

def seperate_bedsummary_by_chromosome(this_species, bed_summary_file):
	"""
	Read a window-ed summary bed file and split it into smaller files according to chromosome
	
	"""
	
	out_filename = [];
	chrom_set = GenomeData.species_chroms[this_species];
	outfile_handler = [0]*len(chrom_set);
	for index in range(0, len(chrom_set)):
		out_filename.append(bed_summary_file+"_temp_"+ chrom_set[index]);
	for index in range(0, len(chrom_set)):
		outfile_handler[index] = open(out_filename[index],'w');
		
	
	infile = open(bed_summary_file,'r');
	for line in infile:					
		""" check to make sure not a header line """
		if not re.match("track", line):
			line_prime = line.strip();
			sline = line_prime.split();
			chrom_index = chrom_set.index(sline[0]);
			(outfile_handler[chrom_index]).write(line);
	for index in range(0, len(chrom_set)):
		(outfile_handler[index]).close();
	
	infile.close();





def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string",
			dest="species", help="species", metavar="<str>")
	parser.add_option("-b", "--bedfile", action="store", type="string",
			dest="bedfile",  metavar="<file>", 
			help="file with summary info in bed graph format")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
		parser.print_help()
		sys.exit(1)

        if GenomeData.species_chroms.has_key(opt.species):
		seperate_bedsummary_by_chromosome(opt.species, opt.bedfile);	
        else:
            	print opt.species + " is not in the species list";




if __name__ == "__main__":
    	main(sys.argv)
