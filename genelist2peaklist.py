#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import Utility
from gene_set_manipulation import *

cat = "/bin/cat"

def main(argv):
	parser = OptionParser()
	parser.add_option("-l", "--genelist", action="store", type="string", dest="genefile", metavar="<file>", help="genelist file")
	parser.add_option("-d", "--dir", action="store", type="string", dest="peakdir", metavar="<file>", help="directory of info")
	parser.add_option("-i", "--infofile", action="store", type="string", dest="infofile", metavar="<file>", help="info file name suffix")
	parser.add_option("-p", "--peakfile", action="store", type="string", dest="peakfile", metavar="<file>", help="peak file")
	parser.add_option("-o", "--outfile", action="store", type="string", dest="outfile", metavar="<file>", help="output file name")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
        	parser.print_help()
        	sys.exit(1)
	
	
	genelist = get_gene_list(opt.genefile, 0)
	l = open("templist.txt", 'w')
	l.close()
	for gene in genelist:
		if Utility.fileExists(opt.peakdir + '/' + gene + opt.infofile):
			file = opt.peakdir + '/' + gene + opt.infofile
			tmpFile = 'templist.txt'
			try:
				if os.system('%s %s >> %s' % (cat, file, tmpFile)): raise
			except: sys.stderr.write(gene + " does not exist\n");
	peaklist = get_gene_list("templist.txt", 0)
	o = open(opt.outfile, 'w')
	f = open(opt.peakfile, 'r')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			id = sline[0]
			if id in peaklist:
				o.write('\t'.join(sline) + '\n')
	o.close()


if __name__ == "__main__":
	main(sys.argv)