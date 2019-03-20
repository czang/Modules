#!/usr/bin/env python
# Copyright (c) 2009 GWU & NHLBI, NIH
# Authors: Chongzhi Zang, Weiqun Peng and Keji Zhao
#
# This software is distributable under the terms of the GNU
# General Public License (GPL) v2, the text of which can be found at
# http://www.gnu.org/copyleft/gpl.html. Installing, importing or otherwise
# using this module constitutes acceptance of the terms of this License.
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
# 
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser


def Coords2BED(inputfile, outputfile):
    infile = open(inputfile, 'r');
    outfile = open(outputfile, 'w');
    for line in infile:
        line = line.strip();
        sline = line.split();
	sline1 = sline[0].split(':');
	sline2 = sline1[2].split('-');
	outfile.write(sline1[1] + '\t' + sline2[0] + '\t' + sline2[1] + '\t0\t0\t+\n');
    infile.close();
    outfile.close();


def main(argv):
    parser = OptionParser()
    parser.add_option("-r", "--eland_file", action="store", type="string",
                      dest="elandfile", help="eland result file", metavar="<file>")
    parser.add_option("-o", "--outfile", action="store", type="string",
                      dest="outfile", help="output file name", metavar="<file>")
    (opt, args) = parser.parse_args(argv)
    if len(argv) < 4:
        parser.print_help()
        sys.exit(1)

    Coords2BED(opt.elandfile, opt.outfile);

if __name__ == "__main__":
    main(sys.argv)
