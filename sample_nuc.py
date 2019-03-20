#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
import random
import shutil
import get_total_tag_counts

def main(argv):
	parser = OptionParser();
	parser.add_option("-f", "--rawtagfile", action="store", type="string",
			  dest="raw_bed_meta_name", help="raw bed file",
			  metavar="<file>");
	parser.add_option("-n", "--desirednumberoftags", action="store", type="int",
			  dest="desired_number_tags", help="desired number of tags",
			  metavar="<int>");
	parser.add_option("-o", "--splicedrawtagfile", action="store", type="string",
			  dest="out_file_name", help="spliced raw bed file",
			  metavar="<file>");
	parser.add_option("-p", "--splicedrawtagfile2", action="store", type="string",
			  dest="out_file_name2", help="spliced raw bed file",
			  metavar="<file>");
	(opt, args) = parser.parse_args(argv);
	if len(argv) < 6:
		parser.print_help()
		sys.exit(1)
	random_sample(opt.desired_number_tags, opt.raw_bed_file, opt.out_file_name, opt.out_file_name2);
	total = get_total_tag_counts.get_total_tag_counts(opt.out_file_name);
	print "The number of tags in " + opt.out_file_name + ' is ' + str(total);


if __name__ == "__main__":
	main(sys.argv)