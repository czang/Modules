#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

"""
* identical tags are counted in
* tags on different strands are counted differently
	-> tag start + fragment size/2
	<- tag start - fragment size/2

IMPORTANT: it is assumed that chromosome are already separated

"""

plus = re.compile("\+");
minus = re.compile("\-");


redundant_threshold = 3;

	
def get_bed_coords(file, fragment_size):
    """
    Read in bed values and remove identical read occurences when
    greater than redundant_threshold
    """

    infile = open(file);
    shift = int(round(fragment_size/2));

    plus_have = {};
    minus_have = {};

    """Initialize dictionaries"""
    for line in infile:
        line = line.strip();
        sline = line.split();
        if plus.match(sline[5]):
            plus_have[atoi(sline[1]) + shift] = 0;
        else:
            if atoi(sline[2]) - shift > 0:
                minus_have[atoi(sline[2]) - shift] = 0;
        minus_have[0] = 0;
    infile.close();
    """ The possibility of maximum beyond the end of chromosome is
    still not taken into account, this will slow things down quite a
    bit, and probably isn't necessary """


    """Read in bed values and store the occurrences at positions"""
    infile = open(file);
    for line in infile:
        line = line.strip();
        sline = line.split();
        if plus.match(sline[5]):
            plus_have[atoi(sline[1]) + shift] += 1
        elif minus.match(sline[5]):
            if atoi(sline[2]) - shift > 0:
                minus_have[atoi(sline[2]) - shift] += 1;
            else:
                minus_have[0] += 1;
    infile.close()

    """Remove the redundancy"""
    taglist = [];
    for t in plus_have.keys():
        add = [t] * plus_have[t];
        taglist.extend(add[0:redundant_threshold]);
    for t in minus_have.keys():
        add = [t] * minus_have[t];
        taglist.extend(add[0:redundant_threshold]);

    """
    sort
    """
    taglist.sort();
    return taglist;

    
def Generate_windows_and_count_tags(taglist, chrom, fragment_size, file):
    """
    taglist: sorted list of positions that includes every tag on a
    chromosome fragment_size: the artificial bin size for binning the
    tags bed_vals: a dictionary keyed by the start of tag_containing
    windows, with value being the tag count in the window.
    
    In this function, the bins are set up using an absolute coordinate
    system.  Namely [0, fragment_size-1),[fragment_size,
    2*fragment_size-1),...
    
    The result writen into the file is guaranteed to be already sorted
    within a chromosome.
    """
   
    outfile = open(file, 'w');
    if(len(taglist)>0):
        current_window_start = (taglist[0]/fragment_size)*fragment_size; 
        tag_count_in_current_window = 1;
	    
        for i in range(0, len(taglist)):
            start = (taglist[i]/fragment_size)*fragment_size;		
            if start == current_window_start: tag_count_in_current_window += 1;
            elif start > current_window_start:
                # All the tags in the previous window have been counted 
                current_window_end = current_window_start + fragment_size -1;
                # write the window to file
                outline = chrom + "\t" + str(current_window_start) + \
                          "\t" + str(current_window_end) + "\t" + \
                          str(tag_count_in_current_window) + "\n";
                outfile.write(outline);
                current_window_start = start;
                tag_count_in_current_window = 1;
            else:
                print 'Something is wrong!!!!!!!';
        current_window_end = current_window_start + fragment_size -1;
        outline = chrom + "\t" + str(current_window_start) + \
                  "\t" + str(current_window_end) + "\t" +  \
                  str(tag_count_in_current_window) + "\n";
        outfile.write(outline);
    outfile.close();
 
 

def main(argv):
    parser = OptionParser()
    parser.add_option("-f", "--tagfile", action="store", type="string",
                      dest="tagfile", help="file with tag coords in bed format",
                      metavar="<file>")
    parser.add_option("-c", "--chrom", action="store", type="string",
                      dest="chrom", help="chromosome name for graph",
                      metavar="<string>")
    parser.add_option("-i", "--fragment_size", action="store", type="int",
                      dest="fragment_size", metavar="<int>",
                      help="fragment size")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output file name",
                      metavar="<file>")
    
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 8:
        parser.print_help()
        sys.exit(1)
        

    tag_list = get_bed_coords(opt.tagfile, opt.fragment_size);
    Generate_windows_and_count_tags(tag_list, opt.chrom,
                                    opt.fragment_size, opt.outfile);


if __name__ == "__main__":
	main(sys.argv)


        
