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
	os.system("mkdir /zzh_gpfs02/yanmingchen/2_sample_test/%s"%prj_name)
	os.system("mkdir /zzh_gpfs02/yanmingchen/2_sample_test/%s/origin"%prj_name)
	barcode_file = open("../barcode_all.txt", "rU")
	barcode_dict = {}
	for line in barcode_file:
		print line
		line =  line.replace("\n", "").split("\t")
		print line
		if len(line) == 2:
			barcode_dict[line[0]] = line[1]
	print barcode_dict
	infile = "%s/%s.assembled_trimed.fastq"%(prj_tree.origin, prj_name)
	print infile
	handle = SeqIO.parse(open(infile, "rU"), "fastq")
	#n = 1
	writer = open("/zzh_gpfs02/yanmingchen/2_sample_test/%s/origin/%s.assembled.fastq"%(prj_name, prj_name), "w")
	for entry in handle:
		seq = str(entry.seq)
		flag = 0
		#print entry.id, [entry.id],str(entry.id)
		#sys.exit(0)
		if "M03098:58:000000000-B3NR2:1:1107:15715:17908"  in  str(entry.id):
			print entry, entry.id
		#sys.exit(0)	
		for (key, value) in barcode_dict.items():
			#print key, value
			if entry.id == "M03098:58:000000000-B3NR2:1:1107:15715:17908" :
				print key, prj_name
			if 	key.split("_")[-1] == "1":
				barcode_name = key.split("_")[-2]
			else:
				barcode_name = key.split("_")[-1]
			if barcode_name not in prj_name:
				if entry.id == "M03098:58:000000000-B3NR2:1:1107:15715:17908" :
					print "yes"
				primer_seq = str(Seq(value).upper())
				
				forword_loc = re.finditer(primer_seq, seq, re.I)
				if entry.id == "M03098:58:000000000-B3NR2:1:1107:15715:17908" :
					print [i for i in forword_loc], primer_seq, seq
				for index, i in enumerate(forword_loc):
					if index == 0 :
						#reads_index += 1
						#seqrecord = SeqRecord(Seq(seq), id = "sample_%s_%s"%(key, reads_index), description = "")
						#SeqIO.write(entry, writer, "fastq")
						#get_id_list.append(entry.id)
						#sys.exit(0)
						flag = 1
						if entry.id == "M03098:58:000000000-B3NR2:1:1107:15715:17908" :
							print flag, forword_loc
				reverse_seq = str(entry.reverse_complement().seq)
				forword_loc = re.finditer(primer_seq, reverse_seq, re.I)
				for index, i in enumerate(forword_loc):
					if index == 0 :
						#reads_index += 1
						#seqrecord = SeqRecord(Seq(seq), id = "sample_%s_%s"%(key, reads_index), description = "")
						#SeqIO.write(entry, writer, "fastq")
						#get_id_list.append(entry.id)
						#sys.exit(0)
						flag = 1
						if entry.id == "M03098:58:000000000-B3NR2:1:1107:15715:17908" :
							print flag

		if flag == 0:
			SeqIO.write(entry, writer, "fastq")
	
			
	
		
		
		

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
