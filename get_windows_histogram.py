#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser

import GenomeData
import Background_island_calculation

Dir = os.getcwd();

def get_windows_histogram(species, summary_graph_file, window_size, output_filename):
	"""
		Build the histogram for the windows according to the tag-count in each window. 
		It uses the summary.graph type file  that includes all the chromosomes. 
		Each line of the file should look like:
			chrom start end tag_count
	
	"""
	totalbp = 0.0;
	value_list = GenomeData.species_chrom_lengths[species].values();
	for values in value_list: totalbp += values;
	totalnumber_windows = round(totalbp/window_size);
	totalnumber_occupied_windows = 0.0;
	#print "Total number of windows is " + str(totalnumber_windows);
		
	current_max = 0;
	windows_hist = [0];

	infile = open(summary_graph_file, 'r');
		
        for line in infile:
        	""" check to make sure not a header line """
		if not re.match("track", line):
			line = line.strip();
			sline = line.split();
			assert (len(sline) == 4);
			tags_in_window = atof(sline[3]);
			tags_in_window = int(tags_in_window)
			if tags_in_window > current_max:
				windows_hist += [0]*(tags_in_window-current_max);
				current_max = tags_in_window;
			windows_hist[tags_in_window] += 1.0;
			totalnumber_occupied_windows += 1.0;
			
	windows_hist[0] = totalnumber_windows - totalnumber_occupied_windows;

	#Get total number of tags
	total_number_tags = 0;
	for i in range (0, len(windows_hist)):
		total_number_tags += i*windows_hist[i];
	print "Total number of tags is " + str(total_number_tags);

	average = (total_number_tags*window_size)*1.0/totalbp; 

	outfile = open(output_filename, 'w');
	outline = "tag_cout"+ "\t" +"Observed_histogram"+  "\t" + "log(h[i]/H)" + "\t" + "Random_expec" + "\n";
	outfile.write(outline);
	for i in range (0, len(windows_hist)):
		if ( windows_hist[i]>0):
			outline = str(i) + "\t" + str(windows_hist[i]) +  "\t" + str(log(windows_hist[i]/totalnumber_windows)) +"\t" + str(Background_island_calculation.poisson(i,average)*totalnumber_windows) + "\n";
			outfile.write(outline);
	outfile.close();
	return windows_hist;
	

    
def main(argv):
	
	parser = OptionParser()
	parser.add_option("-s", "--species", action="store", type="string",
                      dest="species", help="species under consideration", metavar="<str>")
	parser.add_option("-b", "--_summary_graph_file", action="store", type="string",
                      dest="summary_graph_file", help="summary graph file of", metavar="<file>")
	parser.add_option("-w", "--window_size", action="store", type="int",
                      dest="window_size", help="window size", metavar="<int>")
	parser.add_option("-o", "--output_filename", action="store", type="string",
                      dest="output_filename", help="output_filename", metavar="<file>")

	(opt, args) = parser.parse_args(argv)
	if len(argv) < 8:
        	parser.print_help()
        	sys.exit(1)
  	if opt.species not in GenomeData.species_chroms.keys():
		print "The species is not recognized!!";
	else:
		get_windows_histogram(opt.species, opt.summary_graph_file, opt.window_size, opt.output_filename);
    

if __name__ == "__main__":
	main(sys.argv)


        
