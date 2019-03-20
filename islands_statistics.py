#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

import get_total_tag_counts
#from species_chroms import *
import GenomeData


def output_histogram(histogram, filename):
	outfile = open(filename, 'w');
	for i in range(0, len(histogram)):
		if histogram[i]>0:
			outline = str(i) + "\t " + str(histogram[i]) + "\n";
			outfile.write(outline);
	outfile.close();

def find_islands_length_histogram_from_BED(bed_vals, outfilename):
	"""
	For use of obtaining the length histogram of islands. 
	"""
	current_max = 0;
	total_number_islands = 0.0;
	total_length = 0;
	histogram = [0];

	for chrom in bed_val.keys():
		for bed_item in bed_val[chrom]:
			length = bed_item.end - bed_item.start + 1;
			assert (length >= 0);
			total_length += length;
			if length > current_max:
				histogram += [0]*(length-current_max);
				current_max = length;
			histogram[length] += 1.0;
			total_number_islands += 1.0;

	
	outfile = open(outfilename, 'w');
	outline = "# The totoal length of the islands is: " + str(total_length) + "\n";
	outfile.write(outline);
	
	normalization =0.0;
	for i in range(0, len(histogram)):
		if histogram[i]>0:
			normalization += (float(histogram[i])/total_number_islands);
			outline = str(i) + "\t " + str(histogram[i]) + "\n";
			outfile.write(outline);
	outfile.close();
	assert (fabs(normalization-1)<.00000000001);
	return total_length;

"""
Find the island length histogram, 

Return the total length of the islands.
"""
def find_islands_length_histogram_from_file (islands_file, outfilename):
	"""
	For use of obtaining the length histogram of islands. 
	"""
	current_max = 0;
	total_number_islands = 0.0;
	total_length = 0;
	histogram = [0];
	
	infile = open(islands_file, 'r');
	for line in infile:
		if not re.match("track", line):
			line = line.strip();
			sline = line.split();
			assert (len(sline) == 4);
			length = atoi(sline[2]) - atoi(sline[1]) + 1;
			assert (length >= 0);
			total_length += length;
			if length > current_max:
				histogram += [0]*(length-current_max);
				current_max = length;
			histogram[length] += 1.0;
			total_number_islands += 1.0;
	infile.close();
	
	normalization =0.0;
	outfile = open(outfilename, 'w');
	outline = "# The totoal length of the islands is " + str(total_length)+ "\n";
	outfile.write(outline);
	for i in range(0, len(histogram)):
		if histogram[i]>0:
			normalization += (float(histogram[i])/total_number_islands);
			outline = str(i) + "\t " + str(histogram[i]) + "\n";
			outfile.write(outline);
	outfile.close();
	assert (fabs(normalization-1)<.0000000001);
	
	return total_length;

	
def find_islands_histogram_from_bed(bed_val):
	
	current_max = 0;
	total_number_islands = 0.0;
	histogram = [0];
	
	for chrom in bed_val.keys():
		for item in bed_val[chrom]:
			tags_in_island = item.value;
			if tags_in_island > current_max:
				histogram += [0]*int(round((tags_in_island-current_max)));
				current_max = tags_in_island;
			histogram[int(round((tags_in_island)))] += 1.0;
			total_number_islands += 1.0;
	print "Total number of islands is :",  total_number_islands;
	
	return histogram;

	

		
def find_region_histogram (species, total_tag_counts, islands_file, outfilename):
	"""
		Find the histogram of regions given the island file 
		output results to a file following the islands_file convention.
	"""
	
	total_tag_counts_on_islands = get_total_tag_counts.get_total_tag_counts_bed_graph(islands_file);
	genome_length = sum ( (GenomeData.species_chrom_lengths[species]).values() );
	tag_density = float(total_tag_counts)/genome_length;		

	current_max = 0;
	total_number_islands = 0.0;
	histogram = [0];

	infile = open(islands_file, 'r');
	
	for line in infile:
		""" check to make sure not a header line """
		if not re.match("track", line):
			line = line.strip();
			sline = line.split();
			assert (len(sline) == 4);
			tags_in_island = atof(sline[3]);
			if tags_in_island > current_max:
				histogram += [0]*int(round((tags_in_island-current_max)));
				current_max = tags_in_island;
			histogram[int(round((tags_in_island)))] += 1.0;
			total_number_islands += 1.0;
	infile.close();
				
	outfile = open(outfilename, 'w');
	
	outline = "# This is  "+ species  + "\n";
	outfile.write(outline);	
	outline = "# Total number of chromosomes: "+ str(len(GenomeData.species_chrom_lengths[species])) + "\n";
	outfile.write(outline);
	outline = "# Total length of genome is "+ str(genome_length) + "\n";
	outfile.write(outline);
	outline = "# Total number of tags is "+ str(total_tag_counts) + "\n";
	outfile.write(outline);		
	outline = "# Total number of islands is "+ str(total_number_islands) + "\n";
	outfile.write(outline);
	outline = "# Total number of tags on islands is "+ str(total_tag_counts_on_islands) + "\n";
	outfile.write(outline);
	outline = "# Islands tags vs total tags is " + str(total_tag_counts_on_islands*1.0/total_tag_counts) + "\n";
	outfile.write(outline);
	outline = "## of tags" + "\t " + "# of islands" + "\t" + "log((histogram[i])/total_tag_counts)" + "\n";				
	outfile.write(outline);
	
	for i in range(0, len(histogram)):
		if histogram[i]>0:
			outline = str(i) + "\t " + str(histogram[i]) + "\t" + \
				str(log(float(histogram[i])/total_tag_counts)) + "\n";
			outfile.write(outline);
	outfile.close();
	
	return total_number_islands;
	
def main(argv):
	
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="mm8, hg18, background, etc", metavar="<str>")
	parser.add_option("-b", "--summary_graph_file", action="store", type="string",
                      dest="summary_graph_file", help="bed summary graph file of", metavar="<file>")
	parser.add_option("-i", "--islands_file", action="store", type="string",
                      dest="islands_file", help="islands file", metavar="<file>")
	parser.add_option("-o", "--islands_histogram_file", action="store", type="string",
                      dest="islands_histogram_file", help="islands histogram file", metavar="<file>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)

	#tempname = ((opt.islands_file).split('/'))[-1];
	#outfilename = tempname + "_histogram" +".dat"; 

	if opt.species in GenomeData.species_chroms.keys():
		genome_length = sum ( GenomeData.species_chrom_lengths[opt.species].values());
		total_tag_counts = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.summary_graph_file);
		total_number_islands = find_region_histogram(opt.species, total_tag_counts, opt.islands_file, opt.islands_histogram_file);
		
		total_length = find_islands_length_histogram_from_file(opt.islands_file, "length" + opt.islands_histogram_file);
		print "Total islands length is: ", total_length;
		print "Length coverage = total_length_of_islands/genome_length is: ", total_length*1.0/(genome_length - 792451276.0);
		
    	else: 
		print "This species is not in my list!"; 

if __name__ == "__main__":
	main(sys.argv)

