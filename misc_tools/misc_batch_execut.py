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
	for file_name in os.listdir(prj_folder):
		print file_name
		if os.path.isdir(file_name):# and( "K" in file_name or "L" in file_name):
			os.chdir(file_name)
			print os.getcwd()#, os.system("which 2.0.py")
			#os.system("bsub -n 2 -q zzh -e err.6 -o out.6 -P %s \"misc_del_mutibarcodeprimer.py\""%("del" ))
			#os.system("grep \"GCGTAGTACTTACCTGAGGAGACGGTGACC\" ./origin/*.assembled.fastq |wc")
			#os.system("bsub -n 1 -q zzh -e err.0 -o out.0 \"misc-histogram-of-seq-length.py\"")
			#os.system("bsub -n 1 -q zzh -e err.0 -o out.0 \"fastqc %s.R1.fastq %s.R2.fastq\""%(file_name, file_name))
			#os.system("bsub -n 8 -q zzh -e err.7 -o out.7 -P %s \"misc_coverage_distribution.py\""%("cover" ))
			#os.system("bsub -n 8 -q zzh -e err.2 -o out.2 \"2.0.py\"")
			os.system("bsub -n 8 -q zzh -e err.3 -o out.3 -P %s \"misc_clone_distribution.py\""%("clone" ))
			#os.system("bsub -n 8 -q zzh -e err.4 -o out.4 -P %s \"misc_recombanationz_distribution.py\""%("recob" ))
			#os.system("bsub -n 8 -q zzh -e err.5 -o out.5 -P %s \"misc_productive_rate.py\""%("pro" ))
			os.chdir(prj_folder)
		
		

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
