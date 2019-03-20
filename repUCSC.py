# Author :  Shane Bretz
#
# Adapted from USCS.py (Dustin Schones and Keji Zhao)


"""
This module contains the classes to deal with UCSC repetitive
element files.
"""

import re, os, shutil, time, sys, operator
from math import *
from string import *
from optparse import OptionParser

plus = re.compile('\+');
minus = re.compile('\-');

import BED;

repUCSCError = "Error in RepElement class";



class RepElement:
    """
    Class for storing repetitive element information.
    """

    def __init__(self, bin, swScore, milliDiv, milliDel, milliIns, 
                 genoName, genoStart, genoEnd, genoLeft, strand,
                 repName, repClass, repFamily, repStart, repEnd,
                 repLeft, id, age, quality, length):

        self.bin = bin;
        self.swScore = swScore;
        self.milliDiv = milliDiv;
        self.milliDel = milliDel;
        self.milliIns = milliIns;
        self.genoName = genoName;
        self.genoStart = genoStart;
        self.genoEnd = genoEnd;
        self.genoLeft = genoLeft;
        self.strand = strand;
        self.repName = repName;
        self.repClass = repClass;
        self.repFamily = repFamily;
        self.repStart = repStart;
        self.repEnd = repEnd;
        self.repLeft = repLeft;
        self.id = id;
        self.age = age;
        self.quality = quality;
        self.length = length;

    def __setitem__(self, bin, swScore, milliDiv, milliDel, milliIns,
                genoName, genoStart, genoEnd, genoLeft, strand,
                repName, repClass, repFamily, repStart, repEnd,
                repLeft, id, age, quality, length):

        self.bin = bin;
        self.swScore = swScore;
        self.milliDiv = milliDiv;
        self.milliDel = milliDel;
        self.milliIns = milliIns;
        self.genoName = genoName;
        self.genoStart = genoStart;
        self.genoEnd = genoEnd;
        self.genoLeft = genoLeft;
        self.strand = strand;
        self.repName = repName;
        self.repClass = repClass;
        self.repFamily = repFamily;
        self.repStart = repStart;
        self.repEnd = repEnd;
        self.repLeft = repLeft;
        self.id = id;
        self.age = age;
        self.quality = quality;
        self.length = length;
        
    def getAll(self):
        """
        Prints out all important given information.  Information used for
        sorting or filtering not given (quality, age, length of sequence)
        """
        outstring = str(self.bin) + " " + str(self.swScore) + " " + \
                    str(self.milliDiv) + " " + str(self.milliDel) + " " + \
                    str(self.milliIns) + " " + self.genoName + " " + \
                    str(self.genoStart) + " " + str(self.genoEnd) + " " + \
                    str(self.genoLeft) + " " + self.strand + " " + \
                    self.repName + " " + self.repClass + " " + \
                    self.repFamily + " " + str(self.repStart) + " " + \
                    str(self.repEnd) + " " + str(self.repLeft) + " " + \
                    str(self.id) + '\n';
        try:
            return outstring;
        except:
            sys.stderr.write("No UCSC repetive element information for %s\n" % self)
            return ''

#--------------------------
#--------------------------

class KnownRepElements:
    """
    Class to read in UCSC repetitive element files and
    organize information.
    """

    def __init__(self, file=None):
        """
        Reads in a repElement file and builds a dictionary with
        chromosomes as keys and all other information appended in
        a list as values.
        """
        self.repElement_coords = {};  """Dictionary to store all repElements"""
        infile = open(file);
        for line in infile:
            if not re.match("#", line):
                line = line.strip();   """Strip white space off line""" 
                sline = line.split();  """Split line into individual strings (fields)"""
                """ Check to make sure this chromosome is declared in the dictionary.
                    Field 5 is genoName --> the chromosome """
                if sline[5] not in self.repElement_coords.keys():
                    """If not, declare new list"""
                    self.repElement_coords[sline[5]] = []; 
                    """coord is a list in the dictionary repElement_coords 
                        which each instance of repElement is stored with
                        its information."""
                coord = RepElement( sline[0], atoi(sline[1]), atoi(sline[2]),
                                    atoi(sline[3]), atoi(sline[4]), sline[5], 
                                    atoi(sline[6]), atoi(sline[7]), sline[8],
                                    sline[9], sline[10], sline[11], sline[12],
                                    atoi(sline[13]), atoi( sline[14]),
                                    atoi(sline[15]), sline[16],
                                    (atoi(sline[2]) + atoi(sline[3]) + atoi(sline[4])),
                                    0,0);
                """Check to see if the strand is positive or negative, 
                    and calculate quality accordingly."""
                if re.match(plus, coord.strand):
                    coord.quality = (coord.repEnd - coord.repStart)/(float(coord.repEnd - coord.repLeft));
                elif re.match(minus, coord.strand):
                    coord.quality = (coord.repEnd - coord.repLeft)/(float(coord.repEnd -coord.repStart));
                """Calculate length of match"""
                coord.length = coord.genoEnd - coord.genoStart;
                """Add coords to the dictionary under the indicated chromosome"""
                self.repElement_coords[sline[5]].append(coord);

    def FindFamily(self, family):
        """
        Finds all repetitive elements in a certain family and inserts
        them into a dictionary keyed by chromosome.
        """
        Family_coords = {}; """Dict. to  store all elem. of specified family."""
        for chrom in self.repElement_coords.keys():
            """For each chromosome in the dictionary"""
            repFamily_coords = [];
            """New list to store elements of the same family"""
            for element in self.repElement_coords[chrom]:
                """For each read in the chromosome"""
                if str(element.repFamily) == str(family):
                    """Check to see if read matches desired family"""
                    repFamily_coords.append(element);
                    """Add element to new family list"""
            self.repElement_coords[chrom] = repFamily_coords;
            """Replace list in value for the key 'chrom'"""

    def FindName(self, name):
        """
        Finds all repetitive elements carrying a certain name and
        inserts them into a dictionary keyed by chromosome
        """
        for chrom in self.repElement_coords.keys():
            """For each chromosome in the dictionary"""
            repName_coords = [];
            """New list to store elements of the same name"""
            for element in self.repElement_coords[chrom]:
                """For each read in the chromosome"""
                if str(element.repName) == str(name):
                    """Check to see if read matches desired name"""            
                    repName_coords.append(element);
                    """Add element to new name list"""
            self.repElement_coords[chrom] = repName_coords;
            """Replace list in value for the key 'chrom'"""

    def Filter_Elements(self):
        """
        Removes all repetitive element matches which are under
        250 base pairs.  Returns a filtered dictionary.
        """
                  
        for chrom in self.repElement_coords.keys():
            """New dictionary to store a list of filtered elements"""
            filtered_elements = [];
            """For each chromosome in the dictionary"""
            for element in self.repElement_coords[chrom]:
                """For each read in the chromosome"""
                if int(element.length) > 250:
                    filtered_elements.append(element);
                    """Add element to the dictionary if right lenght"""
            self.repElement_coords[chrom] = filtered_elements;
            """Replace the list in [chrom] witht the filtered list"""
    
    def Sort_by_Age(self):
        """
        Sorts the dictionary by age. Age is a value that is the sum
        of milliDiv, milliDel, and milliIns. Sorts youngest to oldest.
        """
        for chrom in self.repElement_coords.keys():
            self.repElement_coords[chrom].sort(key=operator.attrgetter('age'));
        
                   
    def Sort_by_Quality(self):
        """
        Sorts the dictionary by the quality of the repetitive element.
        Quality is determined by a value representing the ratio of
        the match sequence to the whole repetitive element. Sorts
        highest quality to lowest.
        """
        for chrom in self.repElement_coords.keys():
            self.repElement_coords[chrom].sort(key=operator.attrgetter('quality'));
            self.repElement_coords[chrom].reverse();

    def Average_Age(self):
        """
        Finds average age of all elements in the input file
        """
        age_count = 0;
        counter = 0;
        for chrom in self.repElement_coords.keys():
            for element in self.repElement_coords[chrom]:
                age_count = age_count + int(element.age);
                counter = counter + 1;
        avg_age = age_count / float(counter);
        print avg_age;

    def Average_Length(self):
        """
        Finds average length of all elements in the input file
        """
        length_count = 0;
        counter = 0;
        for chrom in self.repElement_coords.keys():
            for element in self.repElement_coords[chrom]:
                length_count = length_count + int(element.length);
                counter = counter + 1;
        avg_length = length_count / float(counter);
        print avg_length;
    
    def keys(self):
        """
        Return a list of keys.
        """
        return self.repElement_coords.keys()

    def __setitem__(self, name, bedcoord):
        """
        Sets a new gene coord (one not in file)
        """
        self.repElement_coords[name] = bedcoord

    def __getitem__(self, name):
        """
        Returns a bed_val indexed by its name or None if no 
        such bed_val exists
        """
        if self.repElement_coords.has_key(name):
            return self.repElement_coords[name]
        else: raise repUCSCError
    
    def __del__(self):
        """
        Delete, delete;
        """
        self.repElement_coords.clear()

    def __contains__(self, item):
        """
        Returns  mapping iterator
        """
        return self.repElement_coords.has_key(item)

    def __iter__(self):
        """
        Returns mapping iterator
        """
        return self.repElement_coords.iterkeys()

    def __len__(self):
        """
        Returns number of gene_coords
        """
        return len(self.gene_coords)

#-----------------------------------------
#-----------------------------------------

def Process_RawData():
    """
    Takes a raw repetive element data file and splits it by
    the repetitive element name, given in the name list.
    FilterElements() and Sort_by_Age() is called to organize
    the data by age of the strand, and filter all elements
    under 250 base pairs.  Creates a new file for each
    repName.
    """

    """List of repElement names to be parsed"""
    Name_List = ["L1PA2", "L1PA3", "L1PA4", "AluY", "L1PA5", "AluSp", "AluSc",
                 "AluSg", "L1PA7", "L1PA6", "AluSx", "L1PA8", "L1PA10", "L1PA11",
                 "L1PA13", "L1PB1", "L1PA12", "L1PA14", "AluJb", "L1PB2", "L1MA1",
                 "L1PA15","L1MA2", "AluJo", "L1MA3", "L1PA16", "L1PB3", "L1MA5",
                 "L1MC", "L1MA8", "L1MB3", "L1MA4", "L1MA6", "L1MA7", "L1MA9",
                 "L1MB8", "L1MC3", "L1ME1", "L1ME2", "MIR", "L2", "L3"];

    for element in Name_List:
        print element;
        for files in os.listdir('/home/shane/data/rep_element_split'):
            os.chdir('/home/shane/data/rep_element_split');
            main_dict = KnownRepElements(files);
            main_dict.FindName(element);
            os.chdir('/home/shane/data/split_by_repname');
            file = open(element, 'a');
            for chrom in main_dict.keys():
                for read in main_dict[chrom]:
                    file.write(read.getAll());                    
            file.close();

    for files in os.listdir('/home/shane/data/split_by_repname'):
        print files;
        os.chdir('/home/shane/data/split_by_repname');
        main_dict = KnownRepElements(files);
        main_dict.Filter_Elements();
        main_dict.Sort_by_Age();
        os.chdir('/home/shane/data/split_by_repname_sorted');
        file = open(files, 'a');
        for chrom in main_dict.keys():
            for read in main_dict[chrom]:
                file.write(read.getAll());
        file.close();

#----------------------------------------
#----------------------------------------

def main(argv):

    """Lines to parse command line arguments from the file repUCSC.sh"""
    parser = OptionParser()
    parser.add_option("-f", "--file", action="store", type="string",
                      dest="datafile", metavar="<file>",
                      help="text file in raw format")    
    parser.add_option("-o","--operation", action="store", type="string",
                      dest="op_type", metavar="<str>",
                      help="name of operation to be performed")
    
    (opt, args) = parser.parse_args(argv)
    
    if len(argv) < 3:
        parser.print_help()
        sys.exit(1)
    
    if (opt.op_type) == "KnownRepElements":
        ## Add functions from KnownRepElements class
        ## above to accomplish task.
        os.chdir('/home/shane/data/sorted_result')
        for files in os.listdir('/home/shane/data/sorted_result'):
            print files;
            main_dict = KnownRepElements(files);
            main_dict.Average_Age();
            main_dict.Average_Length();        

    if (opt.op_type) == "Process_RawData":
        Process_RawData();
       
if __name__ == "__main__":
    main(sys.argv)
