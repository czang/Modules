#!/usr/bin/env python
# Copyright (c) 2008 The George Washington University
# Authors: Weiqun Peng
#
# This software is distributable under the terms of the GNU General
# Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or
# otherwise using this module constitutes acceptance of the terms of
# this License.
#
# Disclaimer
# 
# This software is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# Comments and/or additions are welcome (send e-mail to:
# wpeng@gwu.edu).
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator
from numpy import *

import BED
import GenomeData
import get_total_tag_counts
import find_islands
import Background_simulation_pr
import Background_island_probscore_statistics


def main(argv):
	"""
	Probability scoring with random background model.
	
	"""
	parser = OptionParser()
	
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="mm8, hg18, background, etc", metavar="<str>")
	parser.add_option("-b", "--summarygraph", action="store",type="string", dest="summarygraph", help="summarygraph", metavar="<file>")
	parser.add_option("-w", "--window_size(bp)", action="store", type="int", dest="window_size", help="window_size(in bps)", metavar="<int>")
	parser.add_option("-g", "--gap_size(bp)", action="store", type="int",  dest="gap", help="gap size (in bps)", metavar="<int>")
	parser.add_option("-m", "--score threshold in window", action="store", type="float", dest="window_score_threshold", help="score threshold in window", metavar="<float>")
	parser.add_option("-e", "--evalue ", action="store", type="float", dest="evalue", help="evalue that determines score threshold for islands", metavar="<float>")
	parser.add_option("-o", "--out_sgraph_file", action="store",type="string", dest="out_sgraph_file", help="output for score summary graph", metavar="<file>")
	parser.add_option("-f", "--out_island_file", action="store", type="string", dest="out_island_file", help="output island file name", metavar="<file>")
	
	(opt, args) = parser.parse_args(argv)
	if len(argv) < 16:
        	parser.print_help()
        	sys.exit(1)

	if opt.species in  GenomeData.species_chroms.keys():
		print "Species: ", opt.species;
		print "Window_size: ", opt.window_size;
		print "Gap size: ", opt.gap;
		print "Window score threshold:", opt.window_score_threshold;
		print "E value is:", opt.evalue;
		
		total_read_count = get_total_tag_counts.get_total_tag_counts_bed_graph(opt.summarygraph);
		print "Total read count:", total_read_count
		genome_length = sum (GenomeData.species_chrom_lengths[opt.species].values());
		if opt.species =="hg18":
				genome_length -=  792451276.0;
			# chr1 length
			#genome_length = 195605173.0;
		average = float(total_read_count) * opt.window_size/genome_length; 
		print "Genome Length: ", genome_length;
		print "window average:", average;
		print;
		print "Generate the enriched probscore summary graph file"; 
		#read in the summary graph file
		bed_val = BED.BED(opt.species, opt.summarygraph, "BED_GRAPH");
		#generate the probscore summary graph file, only care about enrichment
		for chrom in bed_val.keys():
			if len(bed_val[chrom])>0:
				for index in xrange(len(bed_val[chrom])):
					read_count = bed_val[chrom][index].value;
					if ( read_count < average):
						score = 0;
					else:
						prob = Background_simulation_pr.poisson(read_count, average);
						if prob <1e-250:
							score = 1000; #outside of the scale, take an arbitrary number.
						else:
							score = -log(prob);
					bed_val[chrom][index].value = score;
					#print chrom, start, read_count, score;
		
		#write the probscore summary graph file
		Background_simulation_pr.output_bedgraph(bed_val, opt.out_sgraph_file);
		
		print "Filter the summary graph to get rid of windows whose scores are less than window_score_threshold";
		#filter the summary graph to get rid of windows whose scores are less than window_score_threshold
		filtered_bed_val = {};
		for chrom in bed_val.keys():
			if len(bed_val[chrom])>0:
				filtered_bed_val [chrom]= [];
				for item in bed_val[chrom]:
					if item.value>=opt.window_score_threshold:
						filtered_bed_val[chrom].append(item);
		
		#Background_simulation_pr.output_bedgraph(filtered_bed_val, opt.out_sgraph_file+".filtered");
		
		print "Determine threshold from random background"; 
		#determine threshold from random background
		bin_size = 0.001;
		#hist_outfile="L" + str(genome_length) + "_W" +str(opt.window_size) + "_G" +str(opt.gap) +  "_s" +str(opt.window_score_threshold) + "_T"+ str(total_read_count) + "_B" + str(bin_size) +"_calculatedprobscoreisland.hist";
		background = Background_island_probscore_statistics.Background_island_probscore_statistics(total_read_count, opt.window_size, opt.gap, opt.window_score_threshold, genome_length, bin_size);
		score_threshold = background.find_island_threshold(opt.evalue); 
		#background.output_distribution(hist_outfile);
		print "The score threshold is: ", score_threshold;
		
		
		print "Make and write islands";
		#make and write islands
		total_number_islands = 0;
		outputfile = open(opt.out_island_file, 'w');
		for chrom in filtered_bed_val.keys():
			if len(filtered_bed_val[chrom])>0:
				islands =find_islands.find_region(filtered_bed_val, chrom, opt.gap, opt.window_size, 10);
				islands =find_islands.find_region_above_threshold(islands, score_threshold);
				total_number_islands += len(islands);
			for i in islands:
				outline = chrom + " " + str(i.start) + " " + str(i.end) + " " + str(i.value) + "\n";	
				outputfile.write(outline);
		outputfile.close();	
		print "total_number_islands:", total_number_islands;
		
	else:
		print "This species is not in my list!"; 

if __name__ == "__main__":
	main(sys.argv)