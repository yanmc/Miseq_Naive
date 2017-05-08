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
	for germline_type in ["V", "J"]:
		V_genes = get_all_gene_name(get_germline_gene(germline_type, chain_type))
		
		infile = "%s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, prj_name)
		writer = csv.writer(open("%s/%s_%s_%s_coverage_distribution.txt"%(prj_tree.data, prj_name, chain_type, germline_type), "w"), delimiter = "\t")
		print "Processing %s " %infile
		handle = open(infile,"rU")
		reader = csv.reader(handle,delimiter = "\t")
		coverage_rate_list = []
		for line in reader:
			assign_result = MyAlignment(line) #Use my_tools's MyAlignment CLASS
			coverage_rate = assign_result.coverage_rate
			if assign_result.assign_type == germline_type and assign_result.sid.split("*")[0] in V_genes:
				coverage_rate_list.append(coverage_rate)
		coverage_rate_list_counter = list_counter(coverage_rate_list)
		X, Y = [], []
		for (rate, number) in coverage_rate_list_counter.items():
			writer.writerow([rate, (float(number)/float(len(coverage_rate_list) ))* 100])
			
			X.append(rate)
			Y.append((float(number)/float(len(coverage_rate_list) )* 100))
		plt.figure(figsize=(9,6))
		plt.bar(X,Y)
		plt.savefig("%s/%s_%s_%s_coverage_distribution.png"%(prj_tree.figure, prj_name, chain_type, germline_type), dpi =300)
		handle.close()
if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	print time.time()
	start = time.time()
	'''
	if "K" in prj_name:
		chain_type = "K"
		main()
	elif "L" in prj_name:
		chain_type = "L"
		main()
	elif "H" in prj_name:
		chain_type = "H"
		main()
		
	else:
		for chain_type in ("H", "K", "L"):
			print "Processing %s chain..."%chain_type
			main()
	'''
	for chain_type in ("H", "K", "L"):
		print "Processing %s chain..."%chain_type
		main()
	
	
	
	end = time.time()
	print prj_name, end-start
	print "Finished"