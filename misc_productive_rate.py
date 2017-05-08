#!/usr/bin/env python
# encoding: utf-8
"""
2.0-process-igblast-output.py 

Created by Mingchen on 2015-05-15.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.
"""

import sys, os, csv, re, glob, copy, subprocess, time, multiprocessing
import traceback, tempfile
import numpy as np
import matplotlib.pyplot as plt
from stackedBarGraph import StackedBarGrapher
import Bio.Alphabet 
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Process, Manager
from bsub import bsub
from mytools import *
from misc_prepare_pbs import *
from misc_get_trimmed_region import *
try:
    import cPickle as pickle
except ImportError:
    import pickle
def main():
	germline_type = "V"
	infile = "%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name)
	print "Processing %s " %infile
	handle = open(infile,"rU")
	reader = csv.reader(handle,delimiter = "\t")
	coverage_rate_list = []
	productiveNum, total  = 0, 0
	for line in reader:
		#print line
		total += 1
		if len(line) > 1 and line[-2] == "Yes":
			productiveNum += 1
	print "The percentage of productive is %s."%(float(productiveNum)/float(total))
	
if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	print time.time()
	start = time.time()
	main()
	
	
	
	end = time.time()
	print prj_name, end-start
	print "Finished"