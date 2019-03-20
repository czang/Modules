#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
#import BED;
#import UCSC;


plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


modifications = ['H2AK5ac','H2AK9ac','H2AK127me1','H2A_Z','H2BK5ac','H2BK5me1','H2BK12ac','H2BK20ac','H2BK120ac','H3K4ac','H3K4me1','H3K4me2','H3K4me3','H3K9ac','H3K9me1','H3K9me2','H3K9me3','H3K14ac','H3K18ac','H3K23ac','H3K23me2','H3K27ac','H3K27me1','H3K27me2','H3K27me3','H3K36ac','H3K36me1','H3K36me3','H3K79me1','H3K79me2','H3K79me3','H3R2me1','H3R2me2','H4K5ac','H4K8ac','H4K12ac','H4K16ac','H4K20me1','H4K20me3','H4K79me2','H4K91ac','H4R3me2']

stemcellmodifications = ['H2AZ','H3K4me1','H3K4me3','H3K9me1','H3K9me3','H3K27me1','H3K27me3','H3K36me3','H4K20me1','PolII']


def gene_modification(datafile):
	file = open(datafile,'r')
	result = {}
	for line in file:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			result[sline[0]] = atof(sline[1])
	return result


def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--inputfiles", action="store", type="string", dest="infile", metavar="<file>", help="input file name")
	parser.add_option("-o", "--output", action="store", type="string", dest="output", metavar="<file>", help="output file name")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	for celltype in ['CD34','CD36']:
		complete_data = {}
		for mod in stemcellmodifications:
			complete_data[mod] = gene_modification(celltype+'-'+mod+opt.infile)
		f = open(celltype+opt.output,'w')
		f.write('#gene')
		for mod in stemcellmodifications:
			f.write('\t' + mod)
		for gene in complete_data['H2AZ'].keys():
			f.write('\n'+gene)
			for mod in stemcellmodifications:
				if complete_data[mod][gene] > 0.1:
					f.write('\t1')
				else:
					f.write('\t0')
		f.close()


if __name__ == "__main__":
	main(sys.argv)