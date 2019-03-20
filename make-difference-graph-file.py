#!/usr/bin/env python
# Copyright (c) 2007 NHLBI, NIH
# Authors: Dustin E Schones, Weiqun Peng, Chongzhi Zang and Keji Zhao
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



""" Read in the tag coords from separate Solexa results and make a
difference map of two libraries

IMPORTANT: it is assumed that chromosome are already separated

"""

plus = re.compile("\+");
minus = re.compile("\-");


def getBedCoords(file, fragment_size):
    """
    *This takes into account the identical tags
    *Tags on different strands are positioned differently
    -> tag start + fragment_size/2
    <- tag start - fragment_size/2
    The stored positions are not the midpoint rathan than the start
    """
    taglist = []
    infile = open(file);
    shift = int(round(fragment_size/2));
    for line in infile:
        """ check to make sure not a header line """
        if not re.match("track", line):
            line = line.strip();
            sline = line.split();
            if plus.match(sline[5]):
                position = atoi(sline[1]) + shift;
                taglist.append(position);
            elif minus.match(sline[5]):
                position = atoi(sline[2]) - shift;
                position = max(0, position);
                taglist.append(position);
    taglist.sort();
    return taglist;



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



def slideWindowAndCountTags(one_tags, two_tags, chrom, window_size, file):
    half_size = float(window_size)/2;
    outfile = open(file, 'a');
    one_tags.sort();
    two_tags.sort();
    one_start = min(one_tags);
    two_start = min(two_tags);
    data_start = min(one_start, two_start);
    one_end = max(one_tags);
    two_end = max(two_tags);
    data_end = max(one_end, two_end);

    data_length = data_end-data_start;
    num_steps = data_length - int(half_size);
    num_steps /= window_size;
    num_steps += 2;  ## psuedo-hack to make sure we get to end of data
    for x in xrange(num_steps):
        midpoint = x*window_size + half_size;
        start = data_start + int(midpoint) - int(half_size);
        end = start + window_size - 1;
        
        one_count = countTagsInWindow(start, end, one_tags);
        two_count = countTagsInWindow(start, end, two_tags);
    
        fold = 0;
        if two_count != 0:
            fold = float(one_count) / float(two_count);
        else:
            fold = float(one_count);

        if fold > 1:
            outline = chrom + "\t" + str(start) + \
                      "\t" + str(end) + "\t" + str(fold) + "\n";
            outfile.write(outline);
    outfile.close();




def main(argv):
    parser = OptionParser()
    parser.add_option("--onetagfile", action="store", type="string",
                      dest="onetagfile",
                      help="first file with tag coords in bed format",
                      metavar="<file>")
    parser.add_option("--twotagfile", action="store", type="string",
                      dest="twotagfile",
                      help="second file with tag coords in bed format",
                      metavar="<file>")
    parser.add_option("-c", "--chrom", action="store", type="string",
                      dest="chrom", help="chromosome name for graph",
                      metavar="<string>")
    parser.add_option("-w", "--window", action="store", type="int",
                      dest="window_size", help="window size to make summary",
                      metavar="<int>")
    parser.add_option("-i", "--fragment_size", action="store", type="int",
                       dest="fragment_size", help="average size of a fragment after CHIP experiment",
                       metavar="<int>")                 
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output file name",
                      metavar="<file>")
    
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 10:
        parser.print_help()
        sys.exit(1)
    
    one_tags = getBedCoords(opt.onetagfile, opt.fragment_size);
    two_tags = getBedCoords(opt.twotagfile, opt.fragment_size);

    slideWindowAndCountTags(one_tags, two_tags, opt.chrom,
                            opt.window_size, opt.outfile);


if __name__ == "__main__":
    main(sys.argv)


        
