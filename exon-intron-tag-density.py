#!/usr/bin/env python
# Copyright (c) 2008 NHLBI, NIH
# Authors: Chongzhi Zang, Dustin Schones, Weiqun Peng and Keji Zhao
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
# zang@gwmail.gwu.edu).


import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

import my_BED;
import UCSC;
import GenomeData


"""
The plan here is to get a profile for each modification across sets of
genes.

Break the 'gene region' up into:

-5kb ... -1 kb ...

...1st 5% of genes, 2nd 5% of gene,  .... last 5% of gene, ...

... +1 kb ... +5 kb

Get the tag densities in each region.
"""

upstream_dist = 5000;
downstream_dist = 5000;
intergenic_window = 1000;
genic_percent_window = 0.05;

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


#hg18_chrom_length_file = "/home/dschones/UCSC/Data/hg18/chrom_lengths.txt"
#dm2_chrom_length_file = "/home/dschones/UCSC/Data/dm2/dm2_chrom_lengths.txt";


chrom_length = {};
def getChromLengths(sp):
    for i in GenomeData.species_chrom_lengths[sp].keys():
	chrom_length[i] = GenomeData.species_chrom_lengths[sp][i]


def binsearch(list, search):    
    """ Classic binary search. If 'search' is in 'list', its index is
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
  """ tally up the tag_starts scores in the range [start, end]"""
  total_tag_starts_score = 0;
  #tag_starts.sort();
  #positions = [(b.position) for b in tag_starts];

  start_ind = binsearch(tag_starts, start);
  end_ind = binsearch(tag_starts, end);
  if type(start_ind) is not int and re.match("no", start_ind):
      val = atoi(start_ind.split(":")[1])
      start_ind = -val-1;
  elif start_ind == -1:
      start_ind = 0;
  if type(end_ind) is not int and re.match("no", end_ind):
      val = atoi(end_ind.split(":")[1]);
      end_ind = -val-1;
      end_ind -= 1 ## if end is not in list, want to go back one
  elif end_ind == -1:
      end_ind = 0;
  for i in range(start_ind, end_ind + 1):
      ## note: range is end inclusive
      total_tag_starts_score += 1;
  return total_tag_starts_score;


def getExonIntronDensities(gene_coords, bed_vals, num_steps):
    tag_starts = [];
    exon_sums = [];
    intron_sums = [];
    exon_sizes = [];
    intron_sizes = [];
    #num_steps = 1.0 / float(genic_percent_window);
    for i in range(int(num_steps)):
        exon_sums.append(0);
        intron_sums.append(0);
        exon_sizes.append(0);
        intron_sizes.append(0);
    for chrom in bed_vals.keys():
        tag_starts = bed_vals[chrom];
        tag_starts.sort();
            
        if chrom in gene_coords.keys() and len(tag_starts)>0:
            for g in gene_coords[chrom]:
		#print g.name
                steps = min(num_steps, g.exonCount-1)
                #g.exonStarts
                #g.exonEnds
                #exon_Start_List = []
                #exon_End_List = []
                #exon_Starts = g.exonStarts.replace(',',' ')
                #exon_Ends = g.exonEnds.replace(',',' ')
                exon_Starts = g.exonStarts.split(',')
                exon_Ends = g.exonEnds.split(',')
                assert len(exon_Starts) == len(exon_Ends)
                if steps > 0:
                    if plus.match(g.strand):
                        for i in range(0, int(steps)):
                            exon_sums[i] += countTagsInWindow(int(exon_Starts[i]), int(exon_Ends[i]), tag_starts);
                            intron_sums[i] += countTagsInWindow(int(exon_Ends[i]), int(exon_Starts[i+1]), tag_starts);
                            exon_sizes[i] += abs(int(exon_Ends[i]) - int(exon_Starts[i]));
                            intron_sizes[i] += abs(int(exon_Starts[i+1]) - int(exon_Ends[i]));
                    elif minus.match(g.strand):
                        for i in range(0, int(steps)):
                            #print exon_Starts
                            exon_sums[i] += countTagsInWindow(int(exon_Starts[-2-i]), int(exon_Ends[-2-i]), tag_starts);
                            intron_sums[i] += countTagsInWindow(int(exon_Ends[-3-i]), int(exon_Starts[-2-i]), tag_starts);
                            exon_sizes[i] += abs(int(exon_Ends[-2-i]) - int(exon_Starts[-2-i]));
                            intron_sizes[i] += abs(int(exon_Starts[-2-i]) - int(exon_Ends[-3-i]));
    return exon_sums, intron_sums, exon_sizes, intron_sizes;




def main(argv):
    parser = OptionParser()
    parser.add_option("-k", "--known_gene_file", action="store", type="string",
                      dest="genefile", help="file with known gene info", metavar="<file>")
    parser.add_option("-b", "--bedfile", action="store", type="string",
                      dest="bedfile", help="file with tags in bed format", metavar="<file>")
    parser.add_option("-n", "--name", action="store", type="string",
                      dest="name", help="name for plotting", metavar="<str>")
    parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="species", metavar="<str>")
    parser.add_option("-c", "--number", action="store", type="int",
                      dest="exonnumber", help="number of exons", metavar="<int>")


    (opt, args) = parser.parse_args(argv)
    if len(argv) < 10:
        parser.print_help()
        sys.exit(1)

    getChromLengths(opt.species);
    gene_coords = UCSC.KnownGenes(opt.genefile);
    bed_vals = my_BED.Starts(opt.species, opt.bedfile);
    num_tags = bed_vals.getNumVals();

    """ print everything out to a file """
    outfilename = '%s-exon-intron-scores' % opt.name;
    outFile = open(outfilename, 'w');


    (exon_counts, intron_counts, exon_seq_sizes, intron_seq_sizes) = getExonIntronDensities(gene_coords, bed_vals, opt.exonnumber);
    for j in range(len(exon_counts)):
        exon_density = float(exon_counts[j]) / float(exon_seq_sizes[j]);
        intron_density = float(intron_counts[j]) / float(intron_seq_sizes[j]);
        exon_density /= float(num_tags);
        intron_density  /= float(num_tags);
        outline = str(j+1) + " " + str(exon_density) + " " + str(intron_density) + "\n";
        outFile.write(outline);
    outFile.close();


if __name__ == "__main__":
    main(sys.argv)


        
