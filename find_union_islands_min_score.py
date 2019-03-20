#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

import BED
from GenomeData import *
import SeparateByChrom
import Utility

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def write(item, out):
	"""
	write one line into outfile. The file openning and closing is handled by outside. 
	item is a BED3 object
	"""
	#chrom, start, end, name, score, strand
	outline = item.chrom + "\t" + str(item.start) + "\t" + str(item.end) + "\t" + str(item.name) + "\t" + str(item.score) + "\n";
	out.write(outline);	


def union_islands_to_file(islandlist, f):
	islandlist.sort(key=operator.attrgetter('start'));
	current = islandlist[0]
	i = 1
	while i < len(islandlist):
		compare = islandlist[i]
		assert current.start <= compare.start
		if compare.start > current.end:
			write(current, f)
			current = compare
			i += 1
		else: 
			current.end = max(current.end, compare.end)
			if compare.score < current.score:
				current.score = compare.score
				current.name = compare.name
			i += 1
	write(current, f)


def main(argv):
	parser = OptionParser()
	parser.add_option("-a", "--islandfile1", action="store", type="string", dest="islandfile1", metavar="<file>", help="BED file with islands info to union")
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, e.g. mm9 or hg18", metavar="<str>")
	parser.add_option("-o", "--outputfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	
	SeparateByChrom.separateByChrom(chroms, opt.islandfile1, '.island1')
	
	for chrom in chroms: 
		f = open(chrom + '.output', 'w')
		bed_vals_1 = BED.BED(opt.species, chrom+'.island1', "BED6", -1)
		if len(bed_vals_1[chrom]) > 0:
			islandlist = bed_vals_1[chrom]
			union_islands_to_file(islandlist, f)
		f.close()

	SeparateByChrom.combineAllGraphFiles(chroms, '.output', opt.outfile);
	SeparateByChrom.cleanup(chroms, '.output')
	SeparateByChrom.cleanup(chroms, '.island1')

if __name__ == "__main__":
	main(sys.argv)       