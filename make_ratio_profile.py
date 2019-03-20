#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;
import UCSC;
from species_chroms import *
import get_total_tag_counts



"""
This module is used to make a UCSC format gene file that can be read by UCSC.py
The input file can be any format gene file given gene name, chromosome, strand, start position and end position in differnt columns.
The input parameters are the 5 colum number for name, chromosome, strand, txStart, txEnd.
The output is the UCSC readable gene file, where other five colums not used are all set zero.

"""


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--inputfile1", action="store", type="string", dest="input_file_1", metavar="<file>", help="naive graph profile 1")
	parser.add_option("-b", "--inputfile2", action="store", type="string", dest="input_file_2", metavar="<file>", help="activated graph profile 2")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="output_file", metavar="<file>", help="output graph profile")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	
	bed_vals_1 = BED.BED(opt.species, opt.input_file_1, "BED_GRAPH", 0) #islands from modification 1
	bed_vals_2 = BED.BED(opt.species, opt.input_file_2, "BED_GRAPH", 0) #islands from modification 2
	
	for chrom in chroms:
		if chrom in bed_vals_1.keys() and chrom in bed_vals_2.keys():
			islands_1 = bed_vals_1[chrom]
			islands_2 = bed_vals_2[chrom]
			f = open(opt.output_file+'_'+chrom,'w')
			for i in range(0, len(islands_1)):
				if i < len(islands_2):
					if islands_1[i].start == islands_2[i].start:
						f.write(str(islands_1[i].start)+'\t'+str(log(islands_2[i].value / islands_1[i].value,2))+'\n')
			f.close()


if __name__ == "__main__":
	main(sys.argv)
