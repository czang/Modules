#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import BED
from GenomeData import *
import separate_by_chrom


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();

def main(argv):
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="species, mm8, hg18", metavar="<str>")
	parser.add_option("-a", "--chipislandfile", action="store", type="string", dest="chipislandbedfile", metavar="<file>", help="chip island file")
	parser.add_option("-b", "--controlislandfile", action="store", type="string", dest="controlislandbedfile", metavar="<file>", help="chip island file")
	parser.add_option("-i", "--chiplibrarysize", action="store", type="int", dest="chipsize", metavar="<int>", help="chip library size")
	parser.add_option("-j", "--controllibrarysize", action="store", type="int", dest="controlsize", metavar="<int>", help="control library size")
	parser.add_option("-t", "--foldchangethreshold", action="store", type="float", dest="foldchangethreshold", metavar="<int>", help="fold change threshold")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="trueislands", metavar="<file>", help="true island file")
	parser.add_option("-p", "--falseislands", action="store", type="string", dest="falseislands", metavar="<file>", help="false island file")
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 14:
        	parser.print_help()
        	sys.exit(1)
	
	if opt.species in species_chroms.keys():
		chroms = species_chroms[opt.species];
	else:
		print "This species is not recognized, exiting";
		sys.exit(1);
	
	chipislands = BED.BED(opt.species, opt.chipislandbedfile, "BED_GRAPH", 0);
	controlislands = BED.BED(opt.species, opt.controlislandbedfile, "BED_GRAPH", 0);
	
	out = open(opt.trueislands, 'w');
	out2 = open(opt.falseislands, 'w')
	for chrom in chroms:
		if chrom in chipislands.keys():
			assert len(chipislands[chrom])==len(controlislands[chrom]);
			for index in xrange(len(chipislands[chrom])):
				island = (chipislands[chrom])[index];
				chipdensity = island.value/float(opt.chipsize);
				controldensity = (controlislands[chrom])[index].value/float(opt.controlsize);
				if controldensity < 0.1/float(opt.controlsize):
					controldensity = 0.1/float(opt.controlsize)
				foldchange = chipdensity/controldensity;
				
				outline = island.chrom +"\t" + str(island.start)+"\t" + str(island.end)+"\t" +str(island.value) +"\t" + str((controlislands[chrom])[index].value)+"\t" +str(foldchange) + "\n";
				
				
				if  foldchange >= opt.foldchangethreshold:
					out.write(outline);
				else:
					out2.write(outline);
	out.close();
	out2.close();
	

if __name__ == "__main__":
	main(sys.argv)