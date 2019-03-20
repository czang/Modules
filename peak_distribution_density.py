#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import bisect

#import BED
import GenomeData
import Utility
import SeparateByChrom

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();

	
def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--bedfile", action="store", type="string", dest="readfile", metavar="<file>", help="raw read file in bed format")
	parser.add_option("-w", "--window", action="store", type="int", dest="window", help="window size", metavar="<int>", default=10000)
	parser.add_option("-t", "--step", action="store", type="int", dest="step", help="sliding step size", metavar="<int>", default=2000)
	parser.add_option("-o", "--outfile", action="store", type="string", dest="out_file", metavar="<file>", help="island read count file")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 6:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in GenomeData.species_chroms.keys():
		chroms = GenomeData.species_chroms[opt.species]
		chrom_lengths = GenomeData.species_chrom_lengths[opt.species]
	else:
		print "This species is not recognized, exiting"
		sys.exit(1);
	
	window = opt.window
	step = opt.step
	
	SeparateByChrom.separateByChrom(chroms, opt.readfile, '.bed1')
	o = open(opt.out_file, 'w')
	for chrom in chroms:
		bedfile = chrom + ".bed1"
		if Utility.fileExists(bedfile):
			print chrom
			position_list = []
			f = open(bedfile, 'r')
			for line in f:
				if not re.match("#", line):
					line = line.strip()
					sline = line.split()
					position = (atoi(sline[1]) + atoi(sline[2])) / 2
					position_list.append(position)
			f.close()
			position_list.sort()
			
			limit = chrom_lengths[chrom]
			start = 1
			end = start + window - 1
			while end < limit: 
				s = bisect.bisect_left(position_list, start)
				e = bisect.bisect_right(position_list, end)
				count = e - s
				o.write(str(start) + '\t' + str(count) + '\n')
				start += step
				end = start + window - 1
	o.close()
	SeparateByChrom.cleanup(chroms, '.bed1')


if __name__ == "__main__":
	main(sys.argv)