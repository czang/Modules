#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
from get_total_tag_counts import *

"""
The majority of the tags are not on a island. 
Example:

mm8-Th-Naive-H3K4me3.bed 	total number of tags is 1326351
		tags on islands ( gap 0, minimum number of tags per window is 2)  is 561162.0
"""


def get_total_tag_counts_bed_graph(bed_graph_file):
	"""
		Get total tag counts given the current experimental run
		file should be a summary bed graph file.
	"""
	counts =0.0;
	infile = open(bed_graph_file,'r');
	for line in infile:
        	""" check to make sure not a header line """
        	if not re.match("track", line):
			line = line.strip();
            		sline = line.split();
			assert ( len(sline) == 4 );	
			counts += atoi(sline[3]);
	infile.close();
	return counts;





def main(argv):
	parser = OptionParser();
	parser.add_option("-f", "--islands_file", action="store", type="string",
			  dest="islands_file", help="file with islands coords in bed format",
			  metavar="<file>");
    
	(opt, args) = parser.parse_args(argv);
	if len(argv) < 2:
		parser.print_help()
		sys.exit(1)

	total = get_total_tag_counts_bed_graph(opt.islands_file);    	
	print "The total number of tags on all the islands is" + opt.islands_file + ' is ' + str(total);
	

if __name__ == "__main__":
	main(sys.argv)