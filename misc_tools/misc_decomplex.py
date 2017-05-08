#!/usr/bin/env python
# encoding: utf-8
"""
3.0.py 

Created by Mingchen on 2015-05-15.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.
"""

import sys, os, csv, re, glob, copy, subprocess, time, multiprocessing
import traceback, tempfile
import pandas as pd
import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
#from stackedBarGraph import StackedBarGrapher
import Bio.Alphabet 
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Process, Manager
from bsub import bsub
from mytools import *
#from misc_prepare_pbs import *
#from misc_get_trimmed_region import *
from collections import Counter
try:
    import cPickle as pickle
except ImportError:
    import pickle
import statsmodels.api as sm
def decpmplex(infiles, primer_seq, key, run_name):
	if os.path.exists("sample_%s"%key):
		pass
	else:
		os.system("mkdir sample_%s"%key)
	get_id_list = []
	for infile in infiles:
		print infile
		handle = SeqIO.parse(open(infile, "rU"), "fastq")
		#n = 1
		for entry in handle:
			#n += 1
			#if n ==1000:
				#print "1000"
				#sys.exit(0)
				#break
			seq = str(entry.seq)
			forword_loc = re.finditer(primer_seq, seq, re.I)
			for index, i in enumerate(forword_loc):
				if index == 0 :
					print i.group()
					sys.exit(0)
					#reads_index += 1
					#seqrecord = SeqRecord(Seq(seq), id = "sample_%s_%s"%(key, reads_index), description = "")
					#SeqIO.write(entry, writer, "fastq")
					get_id_list.append(entry.id)
					#sys.exit(0)
			seq = str(entry.reverse_complement().seq)
			forword_loc = re.finditer(primer_seq, seq, re.I)
			for index, i in enumerate(forword_loc):
				if index == 0 :
					#reads_index += 1
					#seqrecord = SeqRecord(Seq(seq), id = "sample_%s_%s"%(key, reads_index), description = "")
					#SeqIO.write(entry, writer, "fastq")
					get_id_list.append(entry.id)
					#sys.exit(0)
	
	for read_num in ["R1", "R2"]:
		writer = open("./sample_%s/sample_%s.%s.fastq"%(key, key, read_num), "w")
		for infile in infiles:
			if read_num in infile:
				read_num_infile = infile
				print read_num_infile
		read_num_infile_dict = SeqIO.index(read_num_infile, "fastq")
		for entry_id in get_id_list:
			#print entry_id, len(get_id_list)
			SeqIO.write(read_num_infile_dict[entry_id], writer, "fastq")
		#for entry in SeqIO.parse(open(read_num_infile, "rU"), "fastq"):
		
		#	if entry.id in get_id_list:
				#print 
		#		SeqIO.write(entry, writer, "fastq")
				
def main():
	
	#task_pool = Pool(processes = 32)
	for run_name in ["A", "B", "C"]: 
		barcode_file = open("barcode_%s.txt"%run_name, "rU")
		barcode_dict = {}
		for line in barcode_file:
			print line

			line =  line.replace("\n", "").split("\t")
			print line
			if len(line) == 2:
				barcode_dict[line[0]] = line[1]
		print barcode_dict
		#sys.exit(0)
		
		infiles = glob.glob("%s_S*.fastq"%(run_name))
		for (key, value) in barcode_dict.items():
			print key, value
			reads_index = 0
			primer_seq = str(Seq(value).upper())


			#pjobs_ids = task_pool.apply_async(decpmplex, args= (infiles, primer_seq, key, run_name))
			decpmplex(infiles, primer_seq, key, run_name)
	print "Waiting for all subprocesses done..."
	#task_pool.close()
	#task_pool.join()

	print 'All subprocesses done.'
		
		
		

if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	start = time.time()

	main()
	end = time.time()
	print prj_name, end-start
	print "Finished"
