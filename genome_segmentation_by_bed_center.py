#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

from GenomeData import *
import SeparateByChrom
import Utility

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


'''
This module is to segment the genome with a set of boundary regions in BED format, using the center of each BED region as . 

'''


def main(argv):
	parser = OptionParser()
	parser.add_option("-b", "--sample", action="store", type="string", dest="infile", metavar="<file>", help="sample name")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species", metavar="<str>")
	parser.add_option("-o", "--output", action="store", type="string", dest="out_file", metavar="<file>", help="name of output file")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species]
		chrom_lengths = species_chrom_lengths[opt.species]
	else:
		print "This species is not recognized, exiting"
		sys.exit(1)

	if Utility.fileExists(opt.infile):
		SeparateByChrom.separateByChrom(chroms, opt.infile, '.bed1')
	else:
		print opt.bedfile, " not found"
		sys.exit(1)
	
	u = open(opt.out_file, 'w')
	for chrom in chroms:
		if Utility.fileExists(chrom + '.bed1'):
			f = open(chrom + '.bed1', 'r')
			List = []
			for line in f:
				if not re.search("#", line):
					line = line.strip()
					sline = line.split()
					start = atoi(sline[1])
					end = atoi(sline[2])
					center = (start+end) / 2
					List.append(center)
			f.close()
			List.sort()
			length = chrom_lengths[chrom]
			start = 1
			for i in range(0, len(List)):
				end = List[i]
				u.write(chrom + '\t' + str(start) + '\t' + str(end) + '\n')
				start = end + 1
			end = length
			u.write(chrom + '\t' + str(start) + '\t' + str(end) + '\n')
	u.close()
	SeparateByChrom.cleanup(chroms, '.bed1')


if __name__ == "__main__":
	main(sys.argv)