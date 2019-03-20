#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

import gene_set_manipulation

def convert_id(infile, conversion_table, outfile):
	f = open(infile, 'r')
	o = open(outfile, 'w')
	for line in f:
		if not re.match("#", line):
			line = line.strip()
			sline = line.split()
			gene_name = sline[0]
			if gene_name in conversion_table.keys():
				sline[0] = conversion_table[gene_name]	
				outline =  '\t '.join(sline) + '\n'
				o.write(outline)
	f.close()
	o.close()

def main(argv):
    parser = OptionParser()
    parser.add_option("-k", "--conversion_file", action="store", type="string",
                      dest="conversion_file", help="file with conversion table", metavar="<file>")
    parser.add_option("-b", "--gene_file", action="store", type="string",
                      dest="gene_file", help="input file in bed format", metavar="<file>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output file in bed format", metavar="<file>")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 6:
        parser.print_help()
        sys.exit(1)
	
    conversion_table ={}
    conversion_table = gene_set_manipulation.get_conversion_table(opt.conversion_file, 0, 1)
    print len(conversion_table.values()), "elements in the conversion table"
    convert_id(opt.gene_file, conversion_table,  opt.outfile)
    genelist = gene_set_manipulation.get_gene_list(opt.outfile, 0)
    print len(genelist), "genes are converted"
    print len(gene_set_manipulation.find_unique_genes (genelist)), "genes have unique ID"


if __name__ == "__main__":
    main(sys.argv)
