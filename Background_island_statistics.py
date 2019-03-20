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


class Background_island_statistics:
	#  External genomeLength and gapSize are in units of bps
	#  Internal genomeLength and gapSize are in units of windows
	def __init__(self, total_tags, windowSize, gapSize, windowPValue, genomeLength):
		self.tag_density = total_tags * 1.0 / genomeLength;
		self.window_size = windowSize; # In bps. 
		assert(gapSize%windowSize == 0);
		self.gap_size = gapSize/windowSize; # maximum number of windows allowed in a gap 
		self.genome_length = genomeLength/windowSize;
		self.average = self.tag_density * windowSize;
		
		# Precalculate the factorial, poisson, cumulative poisson values up to 250. 
		self.max_index = 250;
		self.fact=[];
		self.poisson_value=[];
		#Cumulative Poisson distribution P_x(<k) = P_x(0)+...+P_x(k-1) 
		self.cumulative_poisson=[];
		self.cumulative_poisson.append(0);
		for index in xrange(self.max_index):
			self.fact.append(self.factorial(index));
			self.poisson_value.append(self.poisson(index, self.average));
			if index > 0:
				self.cumulative_poisson.append(self.cumulative_poisson[index-1] + self.poisson_value[index-1]);
		# get_min_tags_in_window need to use the average, and gap_contribution needs min_tags_in_window
		# So the position of this line is critical.
		self.min_tags_in_window = self.obtain_min_tags_in_window(windowPValue);
		self.gap_contribution = self.gap_factor();
		self.boundary_contribution = self.boundary();
		
		# The previously calculated island_expectation is recorded
		self.island_expectation=[];
		#Precalculate the background island expectation
		for i in xrange(self.min_tags_in_window):
			self.island_expectation.append(0);
		self.island_expectation.append(self.boundary_contribution * self.poisson_value[self.min_tags_in_window] * self.genome_length);



	def factorial(self, m):
		value = 1.0;
		if m != 0:
			while m != 1:
				value = value*m;
				m = m - 1;
		return value;


	# Return the log of a factorial, using Srinivasa Ramanujan's approximation
	def factln(self, m):
		if m<20:  
			value = 1.0;
			if m != 0:
				while m != 1:
					value = value*m;
					m = m - 1;
			return log(value);
		else:
			return m*log(m) -m + log(m*(1+4*m*(1+2*m)))/6.0 + log(pi)/2;


	def poisson(self, i, average):
		if i<20:
			return exp(-average) * average**i / self.factorial(i);
		else:
			exponent = -average + i*log(average) - self.factln(i);
			return exp(exponent);
	
		
	#Returns the thershold value T given the p-value requirement. Namely, P(T)+P(T+1)+ ... +P(infty) <windowPValue. The cumulative_poisson[k] shall be initialized in init. 
	def obtain_min_tags_in_window (self, windowPValue):
		k = 1;
		while ( (1-self.cumulative_poisson[k])>windowPValue):
			k += 1;
		return k;

	"""
		gap is in the unit of windows. In each window in the gap, the
		window could have 0, 1, min_tags_in_windows-1 tags.
		
		say gap = 1, min_tags_in_window= 2, gap_factor = 1 +
		poission(0,a) + poisson(1, a), where 1 represents no gap,
		poisson(0,a) represents a window with 0 tag, 
		poisson(1,a) represents a window with 1 tag,
		
		The gap contribution from each window is not independent
	"""
	def single_gap_factor(self):
		my_gap_factor=0;
		for i in xrange(self.min_tags_in_window): my_gap_factor +=self.poisson_value[i];
		return my_gap_factor;

	def gap_factor(self):
		if self.gap_size == 0: 
			return 1;
		else:
			i = 1;
			gap_contribution = 1; # contribution from no gap
			my_gap_factor = self.single_gap_factor();
			for i in range(1, self.gap_size+1): gap_contribution += pow(my_gap_factor, i);
			return gap_contribution;

	def boundary(self):
		"""
		The condition for boundary is a continuous region of
		unqualified windows longer than gap
		"""
		temp = 0;
		for i in xrange(self.min_tags_in_window): temp += self.poisson_value[i];
		temp = pow(temp, self.gap_size+1); 
		return temp*temp; # start & end 


	#forward method that memorize the calculated results.
	def background_island_expectation (self, score):
		temp=0;
		current_max_score = len(self.island_expectation)-1;
		if score > current_max_score:
			#index is the number of tags in the added window
			for index in range(current_max_score + 1, score+1):
				temp=0;
				for i in range(self.min_tags_in_window, index):
					temp += self.poisson_value[i] * self.gap_contribution * self.island_expectation[index-i];
				# add the contribution from single-window island
 				temp += self.boundary_contribution * self.poisson_value[index] * self.genome_length;
				self.island_expectation.append(temp);
		return self.island_expectation[score];



	def background_island_expectation_wo_single_window_island (self, score):
		if (score < self.min_tags_in_window ): return 0;
		else:
			temp = self.background_island_expectation(score) - self.boundary_contribution * self.poisson_value[score] * self.genome_length;
			return temp;

	def find_island_threshold(self, p_value_threshold):
		"""
		average is the average number of tags in a window:
		opt.tag_density * opt.window_size
		
		This one allows single-window islands.
		Returns the island threshold
		
		p_value threshold currently takes .1
		"""
		threshold = .0000001;
		expectation = [];
		score = self.min_tags_in_window;
		current_expectation=self.background_island_expectation(score);
	
		while current_expectation > threshold:
			score += 1;
			current_expectation=self.background_island_expectation(score);
			

		total = sum (self.island_expectation);
		p_value = [0]*len(self.island_expectation);
		p_value[0] = total;
		for i in range(1, len(self.island_expectation)):	
			p_value[i] = p_value[i-1] - self.island_expectation[i-1];
		for i in range(1, len(p_value)):	
			if p_value[i] < p_value_threshold:
				island_threshold = i;
				break;
		#for i in xrange(len(p_value)):
		#	print i, self.island_expectation[i], p_value[i];	
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
	parser.add_option("-m", "--p_value for a window", action="store", type="float",
			  dest="window_p_value", help="p value to determine min number of tags allowed in a qualified window",
			  metavar="<float>")
	parser.add_option("-l", "--genome_length", action="store", type="float",
			  dest="genome_length", help="genome length in bp",
			  metavar="<float>") 

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 12:
		parser.print_help()
		sys.exit(1)

	#outfile = "Background_t" + str(opt.tag_density)+ "_w"+ str(opt.window_size) + "_g" + str(opt.gap_size) + "_m" + str(opt.min_tags_in_window) + "_calculation.dat";   
	
	#output= open(outfile,'w');
	tag_density = opt.tag_counts/opt.genome_length;
	
	background = Background_island_statistics(opt.tag_counts, opt.window_size, opt.gap_size, opt.window_p_value, opt.genome_length);
	
	print "# Tag count : ", opt.tag_counts;
	print "# Average density of tag on genome is ",  tag_density;
	print "# Genome length: " + str(opt.genome_length) + " or " + str(int(opt.genome_length/opt.window_size)) + " windows";
	print "# Average number of tags in a window: ",  tag_density * opt.window_size;
	print "# Window p-value is ",  opt.window_p_value;
	print "# Minimum num of tags in a qualified window: ", background.min_tags_in_window;
	print "# Gap size: " + str(opt.gap_size) + " bps" ;
	print "# Chosen p value threshold: " + str(opt.p_value_threshold);
	print "# Num tags on island       " + " expected number of islands   " + " pseudo p_value   ";
	#print background.gap_contribution;
	#print background.boundary_contribution;
	#print background.island_expectation;
	
	print background.find_island_threshold(opt.p_value_threshold);
	
if __name__ == "__main__":
    	main(sys.argv)
