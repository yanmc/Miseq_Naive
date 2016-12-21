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
from collections import Counter
try:
    import cPickle as pickle
except ImportError:
    import pickle
import statsmodels.api as sm
import scipy.stats as scist

def main():
	pickle_file = '%s/%s_ploy_gene_position_dict_dump'%(prj_tree.tmp, prj_name)
	ref_seq_id_name = pickle_file.split('_')[-1]
	pickle_file_handle = open(pickle_file, 'rb')
	pickle_tuple = pickle.load(pickle_file_handle)
	ploy_gene_position_dict = pickle_tuple
	pickle_file_handle.close()
	print ploy_gene_position_dict
	for (key, value) in ploy_gene_position_dict.items():
		for position in value:
			if position != 0:
				print key, position
if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	start = time.time()
	pic_type = "_unique"
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
			main()
	end = time.time()
	print prj_name, end-start
	print "Finished"