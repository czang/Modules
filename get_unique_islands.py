#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import BED
from GenomeData import *
from separate_by_chrom import *


def get_unique_bed(bed_vals, chrom, extension):
	end_dic = {}
	value_dic = {}
	for item in bed_vals[chrom]:
		end_dic[item.start] = item.end
		value_dic[item.start] = item.value
	file = open(chrom+extension, 'w')
	start_list = end_dic.keys()
	start_list.sort()
	for start in start_list:
		file.write(chrom + '\t' + str(start) + '\t' + str(end_dic[start]) + '\t' + str(value_dic[start])+'\n')
	file.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species", metavar="<str>")
	parser.add_option("-i", "--bedfile", action="store", type="string", dest="infile", metavar="<file>", help="input island file in BED format")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output unique island file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	separateByChrom(chroms, opt.infile, '.inbed')
	
	for chrom in chroms:
		bed_vals = BED.BED(opt.species, chrom+'.inbed', "BED_GRAPH", 0)
		get_unique_bed(bed_vals, chrom, '.outbed')
	
	combineAllGraphFiles(chroms, '.outbed', opt.outfile)
	cleanup(chroms, '.inbed')
	cleanup(chroms, '.outbed')


if __name__ == "__main__":
	main(sys.argv)