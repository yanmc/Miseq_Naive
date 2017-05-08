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
from collections import Counter
def main():
	infile = "%s/%s_get_recombanation_info.txt"%(prj_tree.igblast_data, prj_name)
	csv_handle = open(infile,"rU")
	reader = csv.reader(csv_handle,delimiter = "\t")
	total, total1,total2, count = 0, 0,0,0
	for index, line in enumerate(reader):
		if len(line) != 1:
			#count += 1
			recomb_result = Recomb(line)
			if recomb_result.v == "N/A" or recomb_result.j == "N/A":
				recomb_result.set_cover_vj("No")
			else:
				recomb_result.set_cover_vj("Yes")
			if recomb_result.productive == "Yes" :
				#:and recomb_result.v.split('*')[0] in germline_gene_list:
				count += 1
			if recomb_result.cover_vj == "Yes" :

				total += 1

			if recomb_result.frame == "In-frame" :

				total1 += 1

			if recomb_result.codon == "No" :

				total2 += 1
	csv_handle.close()
	print count, index, float(count)/index, total, float(total)/index
	print total1, float(total1)/index , total2, float(total2)/index
	#print real_reads_list, count, index+1
	

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
	#infiles = glob.glob("%s/*_real_reads_Variable_region*.fasta"%prj_tree.reads)
	#for infile in infiles:
	#	tanslate_file(infile)
	
	#infiles = glob.glob("%s/*_real_reads_CDR3*.fasta"%prj_tree.reads)
	#for infile in infiles:
	#	tanslate_file(infile)
	
	
	
	end = time.time()
	print prj_name, end-start
	print "Finished"