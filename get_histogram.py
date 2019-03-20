#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


"""
	This module produce a histogram from a file, which is assumed to be of multiple columns separated by tab or space 
	with no missing items. Comment lines which start with # are ignored. 
	column index starts with zero.

"""
def find_histogram (infile_name, column_index, number_of_windows, outfile_name):	

	histogram = [0] * number_of_windows;
	value_list = [];	

	infile = open(infile_name, 'r');
	for line in infile:
		""" check to make sure not a header line """
		if not re.match("#", line):
			line = line.strip();
			sline = line.split();
			assert (len(sline) > column_index);
			if not sline[column_index] == 'None':
				value_list.append(atof(sline[column_index]));
	infile.close();
	print "Total number of data points in this list is " + str(len(value_list));
	

	minimum = min(value_list);
	maximum = max(value_list);
	interval = (1.0000000001 * maximum-minimum)/(number_of_windows);


	x_index=[];
	for index in range(number_of_windows):
		x_index.append(minimum + interval*(index+.5));# Use midpoint to represent the x value.
	
	for value in value_list:
		histogram[int((value-minimum)/interval)] += 1;
	
	outfile = open(outfile_name, 'w');
	outline = "# This a histogram of the "  + str(column_index) + "th column in file" + infile_name + "\n";
	outfile.write(outline);
	outline = "# The minimum and maximum of this list is ("  + str(minimum) +"," + str(maximum) + ")" + "\n";
	outfile.write(outline);
	outline = "# Total number of windows is " + str(number_of_windows) + "\n";
	outfile.write(outline);
	for i in range(0, len(histogram)):
		outline = str(x_index[i]) + "\t " + str(histogram[i]) + "\n";
		outfile.write(outline);
	outfile.close();

def main(argv):
	parser = OptionParser()
	parser.add_option("-i", "--input_file_name", action="store", type="string",
                      dest="infile_name", help="", metavar="<file>")
	parser.add_option("-c", "--column_index", action="store", type="int",
                      dest="column_index", help="the column used for histogram", metavar="<int>")
	parser.add_option("-n", "--number_of_windows", action="store", type="int",
                      dest="number_of_windows", help="number of windows", metavar="<int>")
	parser.add_option("-o", "--ouput_file_name", action="store", type="string",
                      dest="outfile_name", help="output histogram file name", metavar="<file>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)

	find_histogram(opt.infile_name, opt.column_index, opt.number_of_windows, opt.outfile_name);
    

if __name__ == "__main__":
	main(sys.argv)