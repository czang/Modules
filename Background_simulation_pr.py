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
import random
from numpy import *
import bisect

import BED
import GenomeData
import get_total_tag_counts
import get_windows_histogram
import bed_preprocessing
import find_islands
import islands_statistics

"""
============================================================================================
V2.0, 6/22/2008 significant changes happen:
	dropped the requirement that the genome has to have chromosomes of equal length.

	Driver module for a random background: main
"""
def generate_cumulative_hist(hist_file="", histogram=[], outfile=""):
	"""
	Generate cumulative from histogram
	
	The hist_file stores the histogram as double column 
	The histogram is a list of tuples.
	"""
	
	if hist_file != "":
		hist = [];
		f = open(hist_file, "r");
		for line in f:
			line = line.strip();
                        sline = line.split();
			hist.append((atof(sline[0]),atof(sline[1])));
		f.close();
	elif len(histogram)!=0:
		hist = histogram;
	else:
		print "No input, something wrong!";
		
	histy=[];
	cumulative_hist=[];
	for index in xrange(len(hist)):
		histy.append(hist[index][1]);
	for index in xrange(len(hist)):
		partial_sum=sum(histy[index:len(hist)]); # The end is outside of the index
		cumulative_hist.append((hist[index][0], partial_sum));
		
	if outfile != "":
		outf = open(outfile, "w");
		for index in xrange(len(cumulative_hist)):
			outline = str(cumulative_hist[index][0]) + "\t" + str(histy[index])+ "\t" +str(cumulative_hist[index][1]) + "\n";
			outf.write(outline);
		outf.close();
		
	return cumulative_hist;

def generate_random_tags_on_a_chromosome(chromosome_name, chromosome_length, tag_counts, filename, x):
	"""
	Generate random tags for a given chromosome.
	The tags should not be shifted, ie, fragment size is 0. 
	"""
	random.seed(x);
	outfile = open(filename, 'w');
	read_length=25;
	
	for i in xrange(tag_counts):
		position = random.randint(0, chromosome_length);
		outline = chromosome_name + "\t" + str(position) + "\t" + str(position+read_length-1) + "\t"+"U0"+"\t"+"0"+"\t"+"+ \n";
		outfile.write(outline);
	outfile.close();	


def generate_random_tags_on_a_genome(chrom_length, tag_counts, filename, x):
	"""
	Generate random tags for a genome with a number of chromosomes.
	chrom_length is a dictionary of type similar to those in GenomeData
	The tags should not be shifted, ie, fragment size is 0. 
	"""
	random.seed(x);
	read_length=25;
	number_chroms = len(chrom_length);
	chrom_names =  chrom_length.keys();
	chrom_lengths = chrom_length.values(); # will not be in the order it was defined in dictionary.
	genome_length = sum(chrom_lengths);
	chrom_prob = [float(elem)/float(genome_length) for elem in chrom_lengths];
	#chrom_tag_counts = [int(round(tag_counts*elem)) for elem in chrom_prob];
	
	chrom_cumulative_prob=[];
	temp = 0.0;
	for index in xrange(number_chroms):
		temp +=chrom_prob[index];
		chrom_cumulative_prob.append(temp);
		print chrom_cumulative_prob[index];	

	outfile = open(filename, 'w');
	for i in xrange(tag_counts):
		temp = random.random(); #in [0,1)
		chrom_index = bisect.bisect_right(chrom_cumulative_prob, temp);
		position = random.randint(0, chrom_lengths[chrom_index]);
		outline = chrom_names[chrom_index] + "\t" + str(position) + "\t" + str(position+read_length-1) + "\t"+"U0"+"\t"+"0"+"\t"+"+ \n";
		outfile.write(outline);
		#position = random.randint(0, genome_length);
		#Chromosome index starts from 1, till chromosome_number
		#chromosome_index = int(position/chromosome_length) + 1;
		#position = position%chromosome_length; #position within a chromosome
		#outline = "chr" + str(chromosome_index) + "\t" + str(position) + "\t" + str(position+24) + "\t"+"U0"+"\t"+"0"+"\t"+"+ \n";
		#outfile.write(outline);
	outfile.close();	
		


def generate_random_summary (chrom_length, tag_counts, window_size):
	"""
	Directly generate the random summary graph of the entire genome
	summary is a dictionary, it is not exactly the same as the usual summary graph, as 
	it does not record the start and end position, and on the other hand each window  has a value 

	This has to do with how the reads are distributed, so needs to be used for the entire genome at one time.
	"""
	#random.seed(x);
	number_chroms = len(chrom_length);
	chrom_names =  chrom_length.keys(); # will not be in the order it was defined in dictionary.
	length_list = chrom_length.values(); # will not be in the order it was defined in dictionary.
	num_of_windows_in_a_chrom = [int(ceil(float(elem)/window_size)) for elem in length_list];
	genome_length = sum(length_list);
	chrom_prob = [float(elem)/float(genome_length) for elem in length_list];
	
	chrom_cumulative_prob=[];
	temp = 0.0;
	for index in xrange(number_chroms):
		temp +=chrom_prob[index];
		chrom_cumulative_prob.append(temp);
		print "chrom_cumulative_prob[" , index, "]:", chrom_cumulative_prob[index];
	
	summary={};
	for index in xrange(number_chroms):
		chrom_name = chrom_names[index];
		summary[chrom_name] = [0]*num_of_windows_in_a_chrom[index];	
	
	for i in xrange(int(round(tag_counts))):
		temp = random.random(); #in [0,1)
		chrom_index = bisect.bisect_right(chrom_cumulative_prob, temp);
		chrom_name = chrom_names[chrom_index];
		position = random.randint(0, length_list[chrom_index]);
		window_index = int(position/window_size);
		(summary[chrom_name])[window_index] +=1;
	return summary;


def generate_random_summary_in_probscore (chrom_length, tag_counts, window_size):
	"""
	Only focus on enrichment ! so if value<average, the score is zero!!!!

	Directly generate the random summary graph of the entire genome, the score is not in
	read count, but in -log(P).

	summary is a dictionary, it is not exactly the same as the usual summary graph, as 
	it does not record the start and end position, and on the other hand each window  has a value 

	"""
	genome_length = sum(chrom_length.values());
	average = float(tag_counts)*window_size/genome_length;
	
	# This uses the raw method
	# The index of these lists coressponds.
#	chrom_names =  chrom_length.keys(); 
#	length_list = chrom_length.values();
#	summary = generate_random_summary (chrom_length, tag_counts, window_size);
#	for chrom in chrom_names:
# 		for index in xrange(len(summary[chrom])):
# 			value = summary[chrom][index];
# 			if  value >= average:
# 				score = -log(poisson(value, average));
# 			else:
# 				score = 0;
# 			summary[chrom][index] = score;
#	return summary
	# This uses the poisson approach to directly generate the realizations in summary graph.
	summary={};
	total = 0;
	for chrom in chrom_length.keys():
		length = chrom_length[chrom];
		num_windows = int(ceil(float(length)/window_size));
		summary[chrom]=[];
		for index in xrange(num_windows):
			value = random.poisson(average);
			if value >= average:
				score = -log(poisson(value, average));
			else:
				score = 0;	
			(summary[chrom]).append(score);
			total +=value;
	print tag_counts, total;
	return summary;
	
	
def generate_random_islands_histogram(n, species, tag_counts, window_size, gap, window_p_value, idum):
	"""
	n is the number of trials
	chrom_length is a GenomeData type object
	"""
	islands_bed={};
	histogram =[];
	total_histogram=[];
	window_histogram = [];
	
	chrom_length = GenomeData.species_chrom_lengths[species];
	number_chroms = len(chrom_length);
	chrom_names =  chrom_length.keys(); # will not be in the order it was defined in dictionary.
	chrom_lengths = chrom_length.values(); # will not be in the order it was defined in dictionary.
	genome_length = sum(chrom_lengths);
	
	average = tag_counts*window_size/(genome_length*1.0);
	print "average: ", average;
	# the tag count threshold is inclusive
	tag_count_threshold_in_window = find_islands.find_threshold(window_p_value, average);
	print "tag_count_threshold_in_window: ",tag_count_threshold_in_window;
	
	
	max_island_score = 99;
	max_island_count = 5000;
	islands_number_distrib = zeros((max_island_score, max_island_count), int);
	
	
	random.seed(idum);
	
	for i in xrange(n):
		print;
		print;
		print "sample :", i;
		summary = generate_random_summary (chrom_length, tag_counts, window_size);
		FileName = "Randombackground_simulated_L" + str(genome_length) + "_C" + str(number_chroms) + "_T" + str(tag_counts) + "_W" + str(window_size) + "_G" + str(gap) + "_i" + str(idum) +"_R" + str(i);
		
		#Testing block
		#output_bedgraph(translate_summary_to_summarygraph(summary, window_size, 0), FileName + "_summary.graph") ;
		summarygraph = translate_summary_to_summarygraph(summary, window_size, 1);
		window_histogram = bed_preprocessing.combine_histogram(window_histogram, get_windows_histogram.get_windows_histogram_from_bedsummary(species, summarygraph, window_size));
		
		summarygraph = translate_summary_to_summarygraph(summary, window_size, tag_count_threshold_in_window);
		
		
		for chrom in chrom_names:	
			islands_bed={};
			# find_islands.find_region returns a list rather than a dictionary
			islands_bed[chrom] = find_islands.find_region(summarygraph, chrom, gap, window_size);
			#output_bed use append instead of write.
			#output_bed(islands_bed, FileName + "_M" + str(window_p_value) + "_islands.bed");
			histogram = islands_statistics.find_islands_histogram_from_bed(islands_bed);
			
			# get the distributin of number of islands with score s
			for index in xrange(min(len(histogram), max_island_score)):
				if histogram[index]<=max_island_count:
					islands_number_distrib[index][histogram[index]] +=1;
			
			temp = len(histogram) - len(total_histogram);
			if temp>0:
				total_histogram +=[0]*temp;
			for index in xrange(len(histogram)):
				total_histogram[index] += histogram[index];
	
	
	outfile = open("islands_number_distribution.dat", "w");
	for i in xrange(max_island_score):
		outfile.write("\n\n"+ str(i));
		for j in xrange(max_island_count):
			if islands_number_distrib[i][j] > 0:
				outline = str(j) + "\t" + str(islands_number_distrib[i][j]/float(n)) + "\n";
				outfile.write(outline);
	outfile.close();		
	
	print;
	print;
	window_histogram = [elem*1.0/n for elem in window_histogram];
	get_windows_histogram.output_windows_histogram(window_histogram, average, "summarygraph.hist");
	
	total_histogram = [elem*1.0/n for elem in total_histogram];
	return total_histogram;
	
	
def generate_random_probscoreislands_histogram(n, chrom_length, tag_counts, window_size, gap, window_probscore_threshold, bin_size,  idum):
	"""
	n is the number of trials
	chrom_length is a GenomeData type object {chr1:10000000, chr2:13000000, ...}
	gap is in bp
	
	"""
	MAX=100.0;
	bin_number = int(round(MAX/bin_size)) + 1;
	total_histogram=[];
	for index in xrange(bin_number):
		total_histogram.append(0);

	number_chroms = len(chrom_length);
	chrom_names =  chrom_length.keys(); # will not be in the order it was defined in dictionary.
	chrom_lengths = chrom_length.values(); # will not be in the order it was defined in dictionary.
	genome_length = sum(chrom_lengths);
	
	average = tag_counts*window_size/(genome_length*1.0);
	print "average: ", average;
	# the tag count threshold is inclusive
	print "score_threshold_in_window: ",window_probscore_threshold;
	
	random.seed(idum);
	
	for i in xrange(n):
		print i;
		#summary is generated by random placement of tag_counts reads.
		summary = generate_random_summary_in_probscore (chrom_length, tag_counts, window_size);	
		summarygraph = translate_summary_to_summarygraph(summary, window_size, window_probscore_threshold);
		
		for chrom in chrom_names:	
			# find_islands.find_region returns a list rather than a dictionary
			islands_bed= find_islands.find_region(summarygraph, chrom, gap, window_size);
			for item in islands_bed:
				index = int(item.value/bin_size);
				if index >= len(total_histogram):
					total_histogram += [0] *(index-len(total_histogram)+1)
				total_histogram[index] +=1.0;	
		
	total_histogram = [elem*1.0/n for elem in total_histogram];
	bins=[];
	for index in xrange(len(total_histogram)):
		# [0, bin_size) in first bin, [bin_size, 2*bin_size) in second bin, etc
		bins.append(bin_size*index);
	return (total_histogram,bins);		

	
def generate_random_sample_summarygraphs(n, chrom_length, tag_counts, window_size, filename, idum):
	
	for i in xrange(n):
	 	FileNameBuffer = filename + "_n" + str(i) + "_summary.graph"; 
		print FileNameBuffer;
		summary = generate_random_summary (chrom_length,tag_counts, window_size);
		summarygraph = translate_summary_to_summarygraph(summary, window_size, 1);
		output_bed(summarygraph, FileNameBuffer);	
			
def factorial(m):
	value = 1.0;
	if m != 0:
		while m != 1:
			value = value*m;
			m = m - 1;
	return value;


# Return the log of a factorial, using Srinivasa Ramanujan's approximation
def factln(m):
	if m<20:  
		value = 1.0;
		if m != 0:
			while m != 1:
				value = value*m;
				m = m - 1;
		return log(value);
	else:
		return m*log(m) -m + log(m*(1+4*m*(1+2*m)))/6.0 + log(pi)/2;


def poisson(i, average):
	if i<20:
		return exp(-average) * average**i / factorial(i);
	else:
		exponent = -average + i*log(average) - factln(i);
		return exp(exponent);

			
def translate_summary_to_summarygraph(summary, window_size, score_threshold=0):
	"""
	a window is registered only when its score > score_threhold, to save 
	storage and to comply with convention, those windows with zero tags are not 
	saved in the summary bed. 
	"""
	summarygraph = {};
	chrom_names =  summary.keys();
	for chrom in chrom_names:
		bed_list = [];
		for index in xrange(len(summary[chrom])):	
			value = (summary[chrom])[index];
			if value > score_threshold:
				bed = BED.BED_GRAPH(chrom, index*window_size, (index+1)*window_size-1, value);	
				bed_list.append(bed);
		summarygraph[chrom] = bed_list;
	return summarygraph;	
	
	
def output_bedgraph(bed_val, filename):
	outfile = open(filename, 'w');
	for chrom in bed_val.keys():
		for island in bed_val[chrom]:
			outline = island.chrom + "\t" + str(island.start) + "\t" + str(island.end) + "\t" + str(island.value) + "\n";	
                	outfile.write(outline);
	outfile.close();		

def clear_file(filename):
	output=open(filename, 'w');
	output.write("");
	output.close();

		
def main(argv):
	"""
	This main simulates the island statistics of random background 
	"""
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string", dest="species", help="mm8, hg18, background, etc", metavar="<str>")
	parser.add_option("-t", "--total_tag_counts", action="store", type="long", dest="tag_counts", help="total_tag_counts", metavar="<long>")
	parser.add_option("-w", "--window_size(bp)", action="store", type="int", dest="window_size", help="window_size(bp)", metavar="<int>")
	parser.add_option("-g", "--gap_size(bp)", action="store", type="int", dest="gap", help="gap_size(bp)", metavar="<int>")
	parser.add_option("-m", "--p_value for each window", action="store", type="float",
                      dest="window_p_value", help="p value for each window to determine the minimum number of tags required in each window", metavar="<float>")
	parser.add_option("-n", "--number_of_samples", action="store", type="int", dest="n", help="number_of_samples", metavar="<int>")	     
	parser.add_option("-i", "--init_ran_number", action="store", type="int",dest="idum", help="initial seed of random number", metavar="<int>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 14:
        	parser.print_help()
        	sys.exit(1)
		
	if opt.species in GenomeData.species_chroms.keys():
		genome_length = sum( GenomeData.species_chrom_lengths[opt.species].values())
		chromosome_number=len(GenomeData.species_chroms[opt.species]);
		print "genome_length: ", genome_length;
		print "chromosome_number: ", chromosome_number;
		print "tag_counts: ", opt.tag_counts;
		print "window_size: ", opt.window_size;
		print "Gap size: ", opt.gap;
		print "the number of samples :", opt.n;
		print "idum :", opt.idum;

		total_histogram = generate_random_islands_histogram(opt.n, opt.species, opt.tag_counts, opt.window_size, opt.gap, opt.window_p_value, opt.idum);
		FileName = "Randomreads_simulated_S" +opt.species  + "_T" + str(opt.tag_counts) + "_W" + str(opt.window_size) + "_G" + str(opt.gap) + "_i" + str(opt.idum)+ "_M" + str(opt.window_p_value);
		islands_statistics.output_histogram(total_histogram, FileName  + "_islands_total_histogram.dat");
		
    	else: 
		print "This species is not in my list!"; 


if __name__ == "__main__":
	main(sys.argv)
	
