#!/usr/bin/env python
import re, os, sys, shutil
from math import *   
from string import *
from optparse import OptionParser
import operator


def sliding_window (list, window_size):
	landscape=[];
	half = int(round(window_size/2.0));
	index = half;
	temp = 0.0;
	for index in xrange(window_size):
		temp += list[index];
	temp /= window_size;
	landscape.append((start, temp));
	while index < len(list)-half:
		temp  = temp + list[index+half] - list[index-half];
		landscape.append((index, temp));
		index += 1;
	return landscape;
		


