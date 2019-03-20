import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

sys.path.append("/home/zang/Modules")
from gene_set_manipulation import *

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
	
    conversion_table ={};
    conversion_table = get_conversion_table(opt.conversion_file, 0, 1);
    print len(conversion_table.values());
    convert_geneID_in_file(opt.gene_file, 0, conversion_table,  opt.outfile);
    genelist = get_gene_list(opt.outfile, 0)
    print len(genelist);
    print len(find_unique_genes (genelist));


if __name__ == "__main__":
    main(sys.argv)
