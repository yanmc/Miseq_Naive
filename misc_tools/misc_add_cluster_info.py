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
from scipy.stats.stats import pearsonr
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.ticker import NullFormatter
import math

def main():
	Sample_pool1_folder = "/zzh_gpfs02/sunjing/sample_cds_table/39_sample_variant"
	Sample_pool2_folder = "/zzh_gpfs02/sunjing/sample_cds_table/39_sample"
	sample_name_list, clone_dict_list = [], []
	#process 1_sample
	os.chdir(Sample_pool2_folder)
	os.system("mkdir add_info")
	for file_name in os.listdir(Sample_pool2_folder):
		
		print file_name
		print Sample_pool2_folder + "/" + file_name
		print os.path.isfile(file_name)
		if os.path.isfile(file_name):# and file_name == "Sample1":# and( "K" in file_name or "L" in file_name):
			print Sample_pool2_folder + "/" + file_name
			print os.getcwd()#, os.system("which 2.0.py")
			true_filename = file_name[: file_name.rindex(".")]
			print true_filename
			#sys.exit(0)
			infile = "%s"%file_name
			try: 
				infile_variant = "%s/%s.variant_function"%(Sample_pool1_folder, true_filename)
				all_line1_list = []
				get_position, get_lines = [], []
				for line in csv.reader(open(infile_variant, "rU"), delimiter= "\t"):
					#print line[2:4]
					if line[0] in ["splicing", 'exonic;splicing','ncRNA_splicing'] and line[2:4] not in get_position:
						#print line[2:4]
						#print line
						line.insert(0, "line0")
						#print line
						#print line[3:5]
						#sys.exit(0)
						get_lines.append(line)
						get_position.append(line[3:5])
						#print line[3:5], get_position
				os.system("cp %s ./add_info/%s"%(infile, infile))
				handle = csv.writer(open("./add_info/%s"%infile, "a+"), delimiter = "\t")
				print type(get_lines), get_lines[0]
				for get_line in get_lines:
					handle.writerow(get_line)
			except IOError:
				pass
		#sys.exit(0)	
		#print set(all_line1_list)
		#os.chdir(Sample_pool2_folder)
	
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
