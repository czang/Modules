#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator

## get BED module
import BED;

"""
The plan here is to get the tag density for each modifications
in CTCF islands -- this is going to be very similar to the
'gene-tag-density-profile.py' program
"""

##############################
## beware: if you change these,
## you have to change the plotting coords
##############################
upstream_dist = 5000;
downstream_dist = 5000;



outside_island_window = 300;
island_percent_window = 0.05;
###############################

plus = re.compile("\+");
minus = re.compile("\-");
Dir = os.getcwd();


#CTCF_coords_file = "/home/zang/human/hg18-CD4-DNase.graph"


hg18_chrom_length_file = "/home/dschones/UCSC/Data/hg18/chrom_lengths.txt"
dm2_chrom_length_file = "/home/dschones/UCSC/Data/dm2/dm2_chrom_lengths.txt";


chrom_length = {};
def getChromLengths(sp):
    if re.match("dm2", sp):
        chrom_len_file = dm2_chrom_length_file;
    elif re.match("hg18", sp):
        chrom_len_file = hg18_chrom_length_file;
    infile = open(chrom_len_file);
    for line in infile:
        line = line.strip();
        sline = line.split();
        chrom_length[sline[0]] = sline[1];


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





def getUpstreamDensities(CTCF_coords, bed_vals):
    tag_starts = [];
    tag_sums = [];
    seq_sizes = [];

    num_steps = float(upstream_dist) / float(outside_island_window);

    for i in range(int(num_steps)):
        tag_sums.append(0);
        seq_sizes.append(0);
    
    for chrom in CTCF_coords.keys():
	if chrom in bed_vals.keys():
            island_coords = CTCF_coords[chrom];
        
            """ sort the CTCF island coords on this chrom """
            island_coords.sort(key=operator.attrgetter('start'))
                
            tag_starts = bed_vals[chrom];
            tag_starts.sort();      
	    if len(tag_starts) > 0:

            	for i in island_coords:
                    start = i.start - upstream_dist;
                    end = start + outside_island_window;
                    ind = 0;
                    if start >= 0:
                        tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                        seq_sizes[ind] += abs(end - start);
                        while end < i.start:    ## check this <---  < or <=
                            start += outside_island_window;
                            end += outside_island_window;
                            ind += 1;
                            tag_sums[ind] += countTagsInWindow(start, end, tag_starts);               
                            seq_sizes[ind] += abs(end - start);
    return tag_sums, seq_sizes;



def getDownstreamDensities(CTCF_coords, bed_vals):
    tag_starts = [];
    tag_sums = [];
    seq_sizes = [];

    num_steps = float(downstream_dist) / float(outside_island_window);

    for i in range(int(num_steps)):
        tag_sums.append(0);
        seq_sizes.append(0);
        
    for chrom in CTCF_coords.keys():
	if chrom in bed_vals.keys():
            island_coords = CTCF_coords[chrom];

            """ sort the CTCF island coords on this chrom """
            island_coords.sort(key=operator.attrgetter('start'))

            tag_starts = bed_vals[chrom];
            tag_starts.sort();
	    if len(tag_starts) > 0:

            	for i in island_coords:
                    start = i.end;
                    end = start + outside_island_window;
                    ind = 0;
                    if end <= chrom_length[chrom]:
                        tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                        seq_sizes[ind] += abs(end - start);
                        while end < i.end + downstream_dist:
                            start += outside_island_window;
                            end += outside_island_window;
                            ind += 1;
                            tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                            seq_sizes[ind] += abs(end - start);
    return tag_sums, seq_sizes;



def getIslandDensities(CTCF_coords, bed_vals):
    tag_starts = [];
    tag_sums = [];
    seq_sizes = [];
    num_steps = 1.0 / float(island_percent_window);

    for i in range(int(num_steps)):
        tag_sums.append(0);
        seq_sizes.append(0);
    
    for chrom in CTCF_coords.keys():
        island_coords = CTCF_coords[chrom];

        """ sort the CTCF island coords on this chrom """
        island_coords.sort(key=operator.attrgetter('start'))
        
        tag_starts = bed_vals[chrom];
        tag_starts.sort();
        if len(tag_starts) > 0:
            for i in island_coords:
                island_length = abs(i.end - i.start);
                window_size = island_length * island_percent_window;
            
                start = i.start;
                end = start + int(window_size);
                ind = 0;
                tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                seq_sizes[ind] += abs(end - start);
                for i in range(int(num_steps) - 1):
                    start += int(window_size);
                    end += int(window_size);
                    ind += 1;
                    tag_sums[ind] += countTagsInWindow(start, end, tag_starts);
                    seq_sizes[ind] += abs(end - start);
    return tag_sums, seq_sizes;



def main(argv):
    parser = OptionParser()
    parser.add_option("-b", "--bedfile", action="store", type="string",
                      dest="bedfile", help="file with tags in bed format", metavar="<file>")
    parser.add_option("-n", "--name", action="store", type="string",
                      dest="name", help="name for plotting", metavar="<str>")
    parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="species", metavar="<str>")
    parser.add_option("-k", "--coordsfile", action="store", type="string", dest="coords_file", help="file with island coordinates in bed format", metavar="<file>")

    (opt, args) = parser.parse_args(argv)
    if len(argv) < 8:
        parser.print_help()
        sys.exit(1)

    getChromLengths(opt.species);
    
    CTCF_coords = BED.BED(opt.species, opt.coords_file, "BED3")
    bed_vals = BED.Starts(opt.species, opt.bedfile);

    num_tags = bed_vals.getNumVals();

    """ print everything out to a file """
    outfilename = '%s-tag-densities' % opt.name;
    outFile = open(outfilename, 'w');

    out_ind = 0;

    (up_tag_counts, up_seq_sizes) = getUpstreamDensities(CTCF_coords, bed_vals);
    for j in range(len(up_tag_counts)):
        if up_tag_counts[j] > 0:
            tag_density = float(up_tag_counts[j]) / float(up_seq_sizes[j]);
            tag_density /= float(num_tags);
        else:
            tag_density = 0;
        out_ind += 1;
        outline = str(out_ind) + " " + str(tag_density) + "\n";
        outFile.write(outline);

    (island_tag_counts, island_seq_sizes) = getIslandDensities(CTCF_coords, bed_vals);
    for j in range(len(island_tag_counts)):
        if island_tag_counts[j] > 0:
            tag_density = float(island_tag_counts[j]) / float(island_seq_sizes[j]);
            tag_density /= float(num_tags);
        else:
            tag_density = 0;
        out_ind += 1;
        outline = str(out_ind) + " " + str(tag_density) + "\n";
        outFile.write(outline);

    (down_tag_counts, down_seq_sizes) = getDownstreamDensities(CTCF_coords, bed_vals);
    for j in range(len(down_tag_counts)):
        if down_tag_counts[j] > 0:
            tag_density = float(down_tag_counts[j]) / float(down_seq_sizes[j]);
            tag_density /= float(num_tags);
        else:
            tag_density = 0;
        out_ind += 1;
        outline = str(out_ind) + " " + str(tag_density) + "\n";
        outFile.write(outline);
    

if __name__ == "__main__":
    main(sys.argv)


        
