#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

import BED;
import stats;


""" READ IN TAG COORDS FROM SEQUENCING LIBRARY AND ESTIMATE THE
FRAGMENT SIZE

"""

plus = re.compile("\+");
minus = re.compile("\-");


def binsearch(list, search):    
    """ classic binary search -- If 'search' is in 'list', its index is
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



def breakUpStrands(bed_vals):
    plus_starts = [];
    minus_starts = [];
    for b in bed_vals:
        if plus.match(b.strand):
            plus_starts.append(b.start);
        elif minus.match(b.strand):
            minus_starts.append(b.start);
    return (plus_starts, minus_starts)


def getTagIndicesInWindow(tag_starts, start, end):
    """
    Good version - return a list of the tag indices in the range
    """
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
    return prox_start_indices;



def getModeOfFragSizes(bed_vals, small_size, large_size):
    """ make a list of all possible fragment size distances that are
    within small_size to large_size (a reasonable range) ---- then
    calculate the mode of these distances """
    possible_frag_sizes = [];
    for chrom in bed_vals.keys():
        if len(bed_vals[chrom]) > 0:
            (plus_starts, minus_starts) = breakUpStrands(bed_vals[chrom]);
            plus_starts.sort();
            minus_starts.sort();
            for p in plus_starts:
                minus_indices = getTagIndicesInWindow(minus_starts, p + small_size,
                                                      p + large_size);
                for mi in minus_indices:
                    ## cowardly use of abs
                    dist = abs(minus_starts[mi] - p);
                    possible_frag_sizes.append(dist);
    return stats.mode(possible_frag_sizes);
                    
               
class FragSize:
    """
    Class to get fragment size of ChIP-seq library

    Parameters:
    bedfile = self-explanatory
    species = hg18, mm8, dm3, etc.
    small_size = guess for smallest chromatin frag size
    large_size = guess for largest chromatin frag size
    """

    def __init__(self, bedfile=None, species="hg18"):
        if(bedfile):
            self.beds = BED.BED(species, bedfile, "BED2");

    def getMode(self, small_size, large_size):
        self.fragsize = getModeOfFragSizes(self.beds, small_size, large_size);
        return self.fragsize;
    

