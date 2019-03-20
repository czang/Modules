#!/usr/bin/env python
# Copyright (c) 2007 NHLBI, NIH
# Authors: Dustin Schones and Keji Zhao
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

import BED;
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

intergenic_window = 250;
genic_percent_window = 0.05;

plus = re.compile("\+");
minus = re.compile("\-");

Dir = os.getcwd();


hg18_chrom_length_file = "/home/dschones/UCSC/Data/hg18/hg18_chrom_lengths.txt"
dm2_chrom_length_file = "/home/dschones/UCSC/Data/dm2/dm2_chrom_lengths.txt";


chrom_length = {};

def getChromLengths(sp):
	chrom_length = GenomeData.species_chrom_lengths[sp];
#    	if re.match("dm2", sp):
#		chrom_length = dm2_chrom_lengths;
    #if re.match("dm2", sp):
    #    chrom_len_file = dm2_chrom_length_file;
    #elif re.match("hg18", sp):
    #    chrom_len_file = hg18_chrom_length_file;
    #infile = open(chrom_len_file);
    #for line in infile:
    #    line = line.strip();
    #    sline = line.split();
    #    chrom_length[sline[0]] = atoi(sline[1]);


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



def getUpstreamDensities(gene_coords, bed_vals):
    tag_starts = [];
    tag_sums = [];
    seq_sizes = [];
    num_steps = float(upstream_dist) / float(intergenic_window);
    for i in range(int(num_steps)):
        tag_sums.append(0);
        seq_sizes.append(0);

    for chrom in bed_vals.keys():
        tag_starts = [b.start for b in bed_vals[chrom]];
        tag_starts.sort();

        if chrom in gene_coords.keys():
            for g in gene_coords[chrom]:
                if plus.match(g.strand):
                    start = g.txStart - upstream_dist;
                    end = start + intergenic_window;
                    ind = 0;
                    if start >= 0:
                        tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                        seq_sizes[ind] += abs(end - start);
                        for k in range(int(num_steps-1)):
                            start += intergenic_window;
                            end += intergenic_window;
                            ind += 1;
                            tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                            seq_sizes[ind] += abs(end - start);
                elif minus.match(g.strand):
                    end = g.txEnd + upstream_dist;
                    start = end - intergenic_window;
                    ind = 0;
                    if end <= chrom_lengths[chrom]:
                        tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                        seq_sizes[ind] += abs(end - start);
                        for k in range(int(num_steps-1)):
                            end -= intergenic_window;
                            start -= intergenic_window;
                            ind += 1;
                            tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                            seq_sizes[ind] += abs(end - start);
    return tag_sums, seq_sizes;


def getDownstreamDensities(gene_coords, bed_vals):
    tag_starts = [];
    tag_sums = [];
    seq_sizes = [];
    num_steps = float(downstream_dist) / float(intergenic_window);
    for i in range(int(num_steps)):
        tag_sums.append(0);
        seq_sizes.append(0);

    for chrom in bed_vals.keys():
        tag_starts = [b.start for b in bed_vals[chrom]];
        tag_starts.sort();
        if chrom in gene_coords.keys():
            for g in gene_coords[chrom]:
                if plus.match(g.strand):
                    start = g.txEnd;
                    end = start + intergenic_window;
                    ind = 0;
                    if end <= chrom_length[chrom]:
                        tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                        seq_sizes[ind] += abs(end - start);
                        for k in range(int(num_steps-1)):
                            start += intergenic_window;
                            end += intergenic_window;
                            ind += 1;
                            tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                            seq_sizes[ind] += abs(end - start);
                elif minus.match(g.strand):
                    start = g.txStart - downstream_dist;
                    end = start + intergenic_window;
                    ind = 0;
                    if start >= 0:
                        tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                        seq_sizes[ind] += abs(end - start);
                        for k in range(int(num_steps-1)):
                            start += intergenic_window;
                            end += intergenic_window;
                            ind += 1;
                            tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                            seq_sizes[ind] += abs(end - start);
    return tag_sums, seq_sizes;


def getGenicDensities(gene_coords, bed_vals):
    tag_starts = [];
    tag_sums = [];
    seq_sizes = [];
    num_steps = 1.0 / float(genic_percent_window);
    for i in range(int(num_steps)):
        tag_sums.append(0);
        seq_sizes.append(0);
    for chrom in bed_vals.keys():
        tag_starts = [b.start for b in bed_vals[chrom]];
        tag_starts.sort();
        if chrom in gene_coords.keys():
            for g in gene_coords[chrom]:
                gene_length = abs(g.txStart - g.txEnd);
                window_size = gene_length * genic_percent_window;
                if plus.match(g.strand):
                    start = g.txStart;
                    end = start + window_size;
                    ind = 0;
                    tag_sums[ind] += countTagsInWindow(int(start), int(end), tag_starts);
                    seq_sizes[ind] += abs(end - start);
                    for i in range(int(num_steps) - 1):
                        start += window_size;
                        end += window_size;
                        ind += 1;
                        tag_sums[ind] += countTagsInWindow(int(start), int(end), tag_starts);
                        seq_sizes[ind] += abs(end - start);
                elif minus.match(g.strand):
                    end = g.txEnd;
                    start = end - window_size;
                    ind = 0;
                    tag_sums[ind] += countTagsInWindow(int(start), int(end), tag_starts);
                    seq_sizes[ind] += abs(end - start);
                    for i in range(int(num_steps) - 1):
                        start -= window_size;
                        end -= window_size;
                        ind += 1;
                        tag_sums[ind] += countTagsInWindow(int(start), int(end), tag_starts);
                        seq_sizes[ind] += abs(end - start);
    return tag_sums, seq_sizes;




def main(argv):
    parser = OptionParser()
    parser.add_option("-k", "--known_gene_file", action="store", type="string",
                      dest="genefile", help="file with known gene info", metavar="<file>")
    parser.add_option("-b", "--bedfile", action="store", type="string",
                      dest="bedfile", help="file with tags in bed format", metavar="<file>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="outfile name", metavar="<file>")
    parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="species", metavar="<str>")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 8:
        parser.print_help()
        sys.exit(1)

    getChromLengths(opt.species);
    gene_coords = UCSC.KnownGenes(opt.genefile);

    bed_vals = BED.BED(opt.species, opt.bedfile, "BED2");

    num_tags = bed_vals.getNumVals();

    """ print everything out to a file """
    outFile = open(opt.outfile, 'w');

    (up_tag_counts, up_seq_sizes) = getUpstreamDensities(gene_coords, bed_vals);
    out_ind = 0;
    for j in range(len(up_tag_counts)):
        tag_density = float(up_tag_counts[j]) / float(up_seq_sizes[j]);
        tag_density /= float(num_tags);
        out_ind += 1;
        outline = str(out_ind) + " " + str(tag_density) + "\n";
        outFile.write(outline);

    (gene_tag_counts, gene_seq_sizes) = getGenicDensities(gene_coords, bed_vals);
    for j in range(len(gene_tag_counts)):
        tag_density = float(gene_tag_counts[j]) / float(gene_seq_sizes[j]);
        tag_density /= float(num_tags);
        out_ind += 1;
        outline = str(out_ind) + " " + str(tag_density) + "\n";
        outFile.write(outline);
        
    (down_tag_counts, down_seq_sizes) = getDownstreamDensities(gene_coords, bed_vals);
    for j in range(len(down_tag_counts)):
        tag_density = float(down_tag_counts[j]) / float(down_seq_sizes[j]);
        tag_density /= float(num_tags);
        out_ind += 1;
        outline = str(out_ind) + " " + str(tag_density) + "\n";
        outFile.write(outline);

    outFile.close();


if __name__ == "__main__":
    main(sys.argv)


        
