#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


genomefasta_dir = '/cluster2/czang/genomes/humanhg18/rawgenome/'


def GCcontent(seq):
	G = count(seq, 'G')
	C = count(seq, 'C')
	A = count(seq, 'A')
	T = count(seq, 'T')
	total = G + C + A + T
	if total > 0: 
		return round(float(G+C)/float(total)*100.0, 1)
	else:
		return -1

def scan_fasta_GC(fastafile, outfile):
	o = open(outfile, 'w')
	f = open(fastafile, 'r')
	for line in f:
		if  not re.match('>', line):
			line = line.strip()
			GC = GCcontent(line)
			if GC > 0:
				o.write(str(GC) + '\n')
	f.close()
	o.close()


def main(argv):
	parser = OptionParser()
	parser.add_option("-b", "--bedfile", action="store", type="string", dest="bedfile", metavar="<file>", help="regions to be scanned")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 4:
        	parser.print_help()
        	sys.exit(1)
	
	
	scan_fasta_GC(opt.bedfile + '.fa', opt.outfile)
	

if __name__ == "__main__":
	main(sys.argv)