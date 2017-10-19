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

def main():
	barcode_file = open("barcode_1.txt", "rU")
	barcode_dict = {}
	for line in barcode_file:
		print line
		
		line =  line.replace("\n", "").split("\t")
		print line
		if len(line) == 2:
			barcode_dict[line[0]] = line[1]
	print barcode_dict
	#sys.exit(0)
	
	
	infiles = glob.glob("mix1_S6_L001_R1_001.fastq")
	for infile in infiles:
		print infile
		handle = SeqIO.parse(open(infile, "rU"), "fastq")
		for entry in handle:
			for (key, value) in barcode_dict.items():
				#print key, value
				reads_index = 0
				primer_seq = str(Seq(value).upper())
				try:
					os.mkdir("./Run1")
				except:
					pass
				try:
					os.mkdir("./Run1/sample_%s"%key)
				except:
					pass
				writer = open(".Run1/sample_%s/sample_%s.R1.fastq"%(key, key), "wa+")
				seq = str(entry.seq)
				forword_loc = re.finditer(primer_seq, seq, re.I)
				flag = 0
				for index, i in enumerate(forword_loc):
					if index == 0 :
						reads_index += 1
						#seqrecord = SeqRecord(Seq(seq), id = "sample_%s_%s"%(key, reads_index), description = "")
						SeqIO.write(entry, writer, "fastq")
						#sys.exit(0)
						flag = 1
				seq = str(entry.reverse_complement().seq)
				forword_loc = re.finditer(primer_seq, seq, re.I)
				for index, i in enumerate(forword_loc):
					if index == 0 :
						reads_index += 1
						#seqrecord = SeqRecord(Seq(seq), id = "sample_%s_%s"%(key, reads_index), description = "")
						SeqIO.write(entry, writer, "fastq")
						#sys.exit(0)
						flag = 1
				if flag == 1:
					break
		

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
