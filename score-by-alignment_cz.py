#!/usr/bin/env python
# Copyright (c) 2007 NHLBI, NIH
# Authors: Dustin E Schones and Keji Zhao
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
import bisect

import BED;
import UCSC;
from GenomeData import *
import separate_by_chrom

"""
This is a template for the analysis of tag distribution with respect
to a fixed point, such as the TSSs of known genes.
"""

upstream_length = 10000;
downstream_length = 10000;

window_size = 5;
half_size = float(window_size)/2;

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


def fileExists(f):
	try:
		file = open(f)
	except IOError:
		exists = 0
	else:
		exists = 1
	return exists


def binsearch(list, search):    
    """ classic binary search If 'search' is in 'list', its index is
    returned.  If not, - (ind + 1) is returned, where ind is the index
    where 'search' should be in 'list' """
    right = len(list)
    left = 0
    previous_center = -1
    if search < list[0]:
        return -1
    while 1:
        center = (left + right) / 2
        candidate = list[center]
        if search == candidate:
            return center
        if center == previous_center:
            string = "no:" + str(-2 - center);
            return string;
        elif search < candidate:
            right = center
        else:
            left = center
        previous_center = center


def tryExpanding(tag_starts, prox_start_indices):
   """
   This is necessary because the find tags in window function that is
   being behaves such that if there are duplicates on the window
   borders, only one of these is returned.
   """
   
   expanded_start_indices = [];

   """ First try extending at beginning """
   start_ind = prox_start_indices[0];
   test_start = start_ind;
   start = start_ind;
   start_extend = 1;
   while start_extend == 1:
      test_start -= 1;
      if test_start >= 0:
         if tag_starts[test_start] == tag_starts[start_ind]:
            start = test_start;
         else:
            start_extend = 0;
      else:
         start_extend = 0;

   """ Then try extending at end """
   end_ind = prox_start_indices[-1];
   test_end = end_ind;
   end = end_ind;
   end_extend = 1;
   while end_extend == 1:
      test_end += 1;
      if test_end < len(tag_starts):
         if tag_starts[test_end] == tag_starts[end_ind]:
            end = test_end;
         else:
            end_extend = 0;
      else:
         end_extend = 0;

   for i in range(start, end+1):
      expanded_start_indices.append(i);
   return expanded_start_indices;


def countTagsInWindow(start, end, tag_starts):
    """  Good version of function """
    start_ind = binsearch(tag_starts, start)
    end_ind = binsearch(tag_starts, end);
    num_prox_starts = 0;
    add_to_end = 0;
    if type(start_ind) is not int and re.match("no", start_ind):
        val = atoi(start_ind.split(":")[1])
        start_ind = -val-1;
        add_to_end = 1;
    elif start_ind == -1:
        start_ind = 0;
        add_to_end = 1;
    else:
        num_prox_starts += 1;
    if type(end_ind) is not int and re.match("no", end_ind):
        val = atoi(end_ind.split(":")[1]);
        end_ind = -val-1;
    elif end_ind == -1:
        end_ind = 0;
    else:
        num_prox_starts += 1;
    if add_to_end:
        end_ind += 1;
    ## add the number between the indices
    num_prox_starts += (end_ind - start_ind  - 1);
    prox_start_indices = [];
    for n in range(start_ind, start_ind+num_prox_starts):
        prox_start_indices.append(n);
    if len(prox_start_indices) > 0:
        prox_start_indices = tryExpanding(tag_starts, prox_start_indices)
    return len(prox_start_indices);


def breakUpStrands(bed_vals):
    plus_starts = [];
    minus_starts = [];
    for b in bed_vals:
        if plus.match(b.strand):
            plus_starts.append(b.start);
        elif minus.match(b.strand):
            minus_starts.append(b.start);
    print len(plus_starts), len(minus_starts)
    return (plus_starts, minus_starts)



def getScoresNearPoints(coords, bed_vals):
    num_steps = float(upstream_length + downstream_length)/float(window_size);
    plus_scores = [];
    minus_scores = [];
    for i in range(int(num_steps)):
        plus_scores.append(0);
        minus_scores.append(0);
    for chrom in coords.keys():
        genes = coords[chrom];
        if chrom in bed_vals.keys():
            (plus_starts, minus_starts) = breakUpStrands(bed_vals[chrom]);
            plus_starts.sort();
            minus_starts.sort();
            for g in genes:
                if plus.match(g.strand):
                    data_start = g.txStart - upstream_length;
                    score_ind = 0;
                    for x in xrange(int(num_steps)):
                        midpoint = x*window_size + half_size;
                        start = data_start + midpoint - half_size;
                        end = start + window_size - 1;
			plus_count = 0
			minus_count = 0
			if len(plus_starts)>0:
                        	plus_count = countTagsInWindow(int(start), int(end),
                                                       plus_starts);
			if len(minus_starts)>0:
                        	minus_count = countTagsInWindow(int(start), int(end),
                                                        minus_starts)
                        plus_scores[score_ind] += plus_count;
                        minus_scores[score_ind] += minus_count;
                        score_ind += 1;
                elif minus.match(g.strand):
                    data_start = g.txEnd + upstream_length;
                    score_ind = 0;
                    for x in xrange(int(num_steps)):
                        midpoint = x*window_size + half_size;
                        start = data_start - midpoint + half_size;
                        end = start - window_size + 1;
                        ## notice, swapped strands for tags because
                        ## we're working on the crick strand now
			plus_count = 0
			minus_count = 0
			if len(minus_starts)>0:
                        	plus_count = countTagsInWindow(int(end), int(start),
                                                       minus_starts);
			if len(plus_starts)>0:
                        	minus_count = countTagsInWindow(int(end), int(start),
                                                        plus_starts);
                        plus_scores[score_ind] += plus_count;
                        minus_scores[score_ind] += minus_count;
                        score_ind += 1;
    return (plus_scores, minus_scores);


def printOut(plus_scores, minus_scores, outfilename, normalization_factor):
    num_steps = float(upstream_length + downstream_length)/float(window_size);
    score_ind = 0;
    outfile = open(outfilename, 'w');
    data_start = -upstream_length
    for x in xrange(int(num_steps)):
        midpoint = x*window_size + half_size;
        start = data_start + int(midpoint) - int(half_size);
        end = start + window_size - 1;
        plus_norm_score = float(plus_scores[score_ind]) / float(normalization_factor);
        minus_norm_score = float(minus_scores[score_ind]) / float(normalization_factor);

        outline = str(start) + "\t" + str(plus_norm_score) + "\t" + \
                  str(minus_norm_score) + "\n";

        outfile.write(outline);
        score_ind += 1;
    outfile.close();
    


def main(argv):
    parser = OptionParser()
    parser.add_option("-k", "--known_genes_file", action="store", type="string",
                      dest="known_file", help="file with known genes", metavar="<file>")
    parser.add_option("-b", "--bedfile", action="store", type="string",
                      dest="bedfile", help="file with tags in bed format", metavar="<file>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="outfile name", metavar="<file>")
    parser.add_option("-n", "--normalize", action="store_true", dest="norm",
                      help="normalize by tag number")
    parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="species", metavar="<str>")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 8:
        parser.print_help()
        sys.exit(1)
	
    if opt.species in species_chroms.keys():
	chroms = species_chroms[opt.species];
	chrom_lengths = species_chrom_lengths[opt.species];
    else:
	print "This species is not recognized, exiting";
	sys.exit(1);

    coords = UCSC.KnownGenes(opt.known_file);
    separate_by_chrom.separateByChrom(chroms, opt.bedfile, '.bed')
    num_genes = 0
    num_tags = 0
    scores = {}
    for chrom in chroms:
	if chrom in coords.keys():
		num_genes += len(coords[chrom])
		if fileExists(chrom+'.bed'):
			bed_vals = BED.BED(opt.species, chrom+'.bed', "BED2");
			num_tags += bed_vals.getNumVals()
			scores[chrom] = getScoresNearPoints(coords, bed_vals);
    separate_by_chrom.cleanup(chroms, '.bed')
    num_steps = float(upstream_length + downstream_length)/float(window_size)
    plus_score_profile = [0]*int(num_steps)
    minus_score_profile = [0]*int(num_steps)
    for chrom in scores.keys():
	(plus_scores, minus_scores) = scores[chrom]
	for i in range(0, len(plus_scores)):
		plus_score_profile[i] += plus_scores[i]
		minus_score_profile[i] += minus_scores[i]
    normalization_factor = 1;
    if opt.norm:
        normalization_factor *= num_tags;
	normalization_factor *= num_genes;
	normalization_factor *= window_size;
    
    printOut(plus_score_profile, minus_score_profile,
             opt.outfile, normalization_factor);
    

if __name__ == "__main__":
    main(sys.argv)
