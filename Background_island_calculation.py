#!/usr/bin/python
# Copyright (c) 2007 NHLBI, NIH
# Authors: Weiqun Peng, Chongzhi Zang and Keji Zhao
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
# schonesde@mail.nih.gov).

import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

"""
Program for calculation of background probability for islands, the
method is simple.  Let's use the simplest case of no gap,
minimum-number-of tag in a window for explanation.  An island has to
start with a window without tags, then a continuous stretch of windows
with tags, and end with another window without tags.  The start
position can be anywhere along the genome. The stretch of filled
windows are independent from each other. So the probability (per bp)
of a island of score 2

P(0, a) * [P(1,a)^2 + P(2,a)]* P(0,a)

The number of tags in a window follows a poisson distribution.


This module is for calculating the background probability
distribution.  Background; given tag density, given genome length,
randomly place the tags Probability distribution P(k): the expected
number of island with given score k on a genome.

Weiqun Peng

*The analytical recursive relationship is developed by Chen Zeng and
 Weiqun Peng

# Currently single window with high counts of tags does not count as
# an island!
"""

	
def factorial(m):
	value = 1.0;
	if m != 0:
		while m != 1:
			value = value*m;
			m = m - 1;
	return value;


# Return the log of a factorial
def factln(m):
	if m<20:  
		return log(factorial(m));
	else:
		return m*log(m) -m + log(m*(1+4*m*(1+2*m)))/6.0 + log(pi)/2;
	

def poisson(i, average):
	exponent = -average + i*log(average) - factln(i);
	return exp(exponent);


def gap_factor (gap, min_tags_in_window, average):
	"""
	gap is in the unit of windows. In each window in the gap, the
	window could have 0, 1, min_tags_in_windows-1 tags.
	
	say gap = 1, min_tags_in_window= 2, gap_factor = 1 +
	poission(0,a) + poisson(1, a), where 1 represents no gap,
	poisson(0,a) represents a window with 0 tag, 
	poisson(1,a) represents a window with 1 tag,
	
	The gap contribution from each window is not independent
	"""
	
	assert min_tags_in_window >= 1;
	if gap == 0: return 1;
	else:
		i = 1;
		gap_contribution = 1; # contribution from no gap
		my_gap_factor = 0;
		for i in range(0, min_tags_in_window): my_gap_factor +=poisson(i, average);
		for i in range(1, gap+1): gap_contribution += pow(my_gap_factor, i);
		return gap_contribution;


def boundary(gap, min_tags_in_window, average):
	"""
	The condition for boundary is a continuous region of
	unqualified windows longer than gap
	"""
	assert min_tags_in_window >= 1;
	temp = 0;
	for i in range(0, min_tags_in_window): temp += poisson(i, average);
	temp = pow(temp, gap+1); 
	return temp*temp; # start & end 


def background_island_expectation (score, min_tags_in_window,  gap, average, genome_length):
	"""
	Recursive
	"""
	if (score < min_tags_in_window ): return 0;
	elif (score == min_tags_in_window ) : return boundary(gap, min_tags_in_window, average)* poisson(score, average) * genome_length;
	else:
		temp = 0;
		for index in range(min_tags_in_window, score):
			temp += poisson(index, average) * \
				gap_factor (gap, min_tags_in_window, average) * \
				background_island_expectation(score-index, min_tags_in_window, gap, average, genome_length);
		temp += boundary(gap, min_tags_in_window, average) * poisson(score, average) * genome_length;
		return temp;

def background_island_expectation_wo_single_window_island (score, min_tags_in_window,  gap, average, genome_length):
	if (score < min_tags_in_window ): return 0;
	else:
		temp = background_island_expectation(score, min_tags_in_window,  gap, average, genome_length) - \
		       boundary(gap, min_tags_in_window, average) * poisson(score, average) * genome_length;
		return temp;

def find_island_threshold(p_value_threshold, min_tags_in_window, gap, average, scaled_genome_length):
	"""
	average is the average number of tags in a window:
	opt.tag_density * opt.window_size
	
	This one allows single-window islands.
	Returns the island threshold
	
	p_value threshold currently takes .1
	"""
	threshold = .00001;
	expectation = [];
	score = min_tags_in_window;
	for i in xrange(score): expectation.append(0);
	current_expectation=background_island_expectation(score,
							  min_tags_in_window,
							  gap,
							  average,
							  scaled_genome_length);
	expectation.append(current_expectation);
	
	while current_expectation > threshold:
		score += 1;
		current_expectation=background_island_expectation(score,
								  min_tags_in_window,
								  gap,
								  average,
								  scaled_genome_length);
		expectation.append(current_expectation);
		#print score, current_expectation;

	total = 0.0
	for i in xrange(len(expectation)):
		total += expectation[i];
	p_value = [0]*len(expectation);
	p_value[0] = total;
	for i in range(1, len(expectation)):	
		p_value[i] = p_value[i-1] - expectation[i-1];
	for i in range(1, len(p_value)):	
		if p_value[i] < p_value_threshold:
			island_threshold = i;
			break;
	#for i in xrange(len(p_value)):
		#print i, expectation[i], p_value[i];	
	return island_threshold;


def main(argv):
	parser = OptionParser();
	parser.add_option("-p", "--p_value_threshold", action="store", type="float",
			  dest="p_value_threshold", help="p_value_threshold",
			  metavar="<float>") 
	parser.add_option("-t", "--tag_counts", action="store", type="float",
			  dest="tag_counts", help="tag counts from experimental data",
			  metavar="<float>") 
	parser.add_option("-w", "--window_size", action="store", type="int",
			  dest="window_size", help="window size in bp to make summary",
			  metavar="<int>")
	parser.add_option("-g", "--gap_size", action="store", type="int",
			  dest="gap_size", help="number of gaps allowed in the island in terms of windows",
			  metavar="<int>")
	parser.add_option("-m", "--min_tags_in_window", action="store", type="int",
			  dest="min_tags_in_window", help="min number of tags allowed in a qualified window",
			  metavar="<int>")
	parser.add_option("-l", "--genome_length", action="store", type="float",
			  dest="genome_length", help="genome length in bp",
			  metavar="<float>") 

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 10:
		parser.print_help()
		sys.exit(1)

	#outfile = "Background_t" + str(opt.tag_density)+ "_w"+ str(opt.window_size) + "_g" + str(opt.gap_size) + "_m" + str(opt.min_tags_in_window) + "_calculation.dat";   
	
	#output= open(outfile,'w');
	tag_density = opt.tag_counts/opt.genome_length;
	print "# Average density of tag on genome is ",  tag_density;
	print "# Genome length: " + str(opt.genome_length) + " or " + str(int(opt.genome_length/opt.window_size)) + " windows";
	print "# Minimum num of tags in a qualified window: " + str(opt.min_tags_in_window);
	print "# Gap size: " + str(opt.gap_size) + " windows or ", opt.gap_size*opt.window_size, " bps" ;
	print "# Chosen p value threshold: " + str(opt.p_value_threshold);
	print "# Num tags on island       " + " expected number of islands   " + " pseudo p_value   ";
	
	average = tag_density * opt.window_size;
	scaled_genome_length = int(opt.genome_length/opt.window_size);
	
	print find_island_threshold(opt.p_value_threshold, opt.min_tags_in_window, opt.gap_size, average, scaled_genome_length);
	
if __name__ == "__main__":
    	main(sys.argv)
