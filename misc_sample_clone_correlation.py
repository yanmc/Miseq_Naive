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
def scatter_hist2(prj_folder, sample_name_X, sample_name_Y, x,y):
	nullfmt   = NullFormatter()         # no labels

	# definitions for the axes
	left, width = 0.1, 0.8
	bottom, height = 0.1, 0.8
	bottom_h = left_h = left+width+0.02

	rect_scatter = [left, bottom, width, height]
	#rect_histx = [left, bottom_h, width, 0.2]
	#rect_histy = [left_h, bottom, 0.2, height]

	# start with a rectangular Figure
	fig = plt.figure(1, figsize=(8,8))

	axScatter = plt.axes(rect_scatter)
	#axHistx = plt.axes(rect_histx)
	#axHisty = plt.axes(rect_histy)

	# no labels
	#axHistx.xaxis.set_major_formatter(nullfmt)
	#axHisty.yaxis.set_major_formatter(nullfmt)

	# the scatter plot:
	axScatter.scatter(x, y)

	# now determine nice limits by hand:
	binwidth = 0.25
	xymax = np.max( [np.max(np.fabs(x)), np.max(np.fabs(y))] )
	#lim = ( int(xymax) +1) #* binwidth
	lim =xymax #* binwidth

	axScatter.set_xlim( (0, lim) )
	axScatter.set_ylim( (0, lim) )

	#bins = np.arange(-lim, lim + binwidth, binwidth)
	#axHistx.hist(x, bins=bins)
	#axHisty.hist(y, bins=bins, orientation='horizontal')

	#axHistx.set_xlim( axScatter.get_xlim() )
	#axHisty.set_ylim( axScatter.get_ylim() )

	#plt.show()
	#print help(plt.savefig)
	plt.savefig('%s/%s_%s_clone_correlation.png'%(prj_folder, sample_name_X, sample_name_Y), dpi = 300)
	del fig
	plt.close()
def main():
	Sample_pool1_folder = "/zzh_gpfs02/yanmingchen/2_sample_test"
	Sample_pool2_folder = "/zzh_gpfs02/yanmingchen/1_sample"
	sample_name_list, clone_dict_list = [], []
	#process 1_sample
	os.chdir(Sample_pool2_folder)
	for file_name in os.listdir(Sample_pool2_folder):
		clone_dict = {}
		print file_name
		print Sample_pool2_folder + "/" + file_name
		print os.path.isdir(file_name)
		if os.path.isdir(file_name):# and file_name == "Sample1":# and( "K" in file_name or "L" in file_name):
			print Sample_pool2_folder + "/" + file_name
			os.chdir(Sample_pool2_folder + "/" + file_name)
			print os.getcwd()#, os.system("which 2.0.py")
			prj_tree = ProjectFolders(Sample_pool2_folder + "/" + file_name)
			prj_name = fullpath2last_folder(prj_tree.home)
			try :
				clone_distribution_file = "%s/%s_H_clone_frequency_all.txt"%(prj_tree.data, prj_name)
				
				for line in csv.reader(open(clone_distribution_file, "rU"), delimiter= "\t"):
					if "Frequency" not in line[-1]:
						clone_dict[tuple(line[:3])] = line[4]
				sample_name_list.append(file_name)
				clone_dict_list.append(clone_dict)
			except IOError:
				pass
			
			
			
			os.chdir(Sample_pool2_folder)
	#process 2_sample
	os.chdir(Sample_pool1_folder)
	for file_name in os.listdir(Sample_pool1_folder):
		clone_dict = {}
		print file_name
		print Sample_pool1_folder + "/" + file_name
		print os.path.isdir(file_name)
		if os.path.isdir(file_name) :# and file_name == "Sample1":# and( "K" in file_name or "L" in file_name):
			print Sample_pool1_folder + "/" + file_name
			os.chdir(Sample_pool1_folder + "/" + file_name)
			print os.getcwd()#, os.system("which 2.0.py")
			prj_tree = ProjectFolders(Sample_pool1_folder + "/" + file_name)
			prj_name = fullpath2last_folder(prj_tree.home)
			
			clone_distribution_file = "%s/%s_H_clone_frequency_all2.txt"%(prj_tree.data, prj_name)
			
			for line in csv.reader(open(clone_distribution_file, "rU"), delimiter= "\t"):
				if "Frequency" not in line[-1]:	
					clone_dict[tuple(line[:3])] = line[4]
			sample_name_list.append(file_name)
			clone_dict_list.append(clone_dict)
			
			
			os.chdir(Sample_pool1_folder)
	print sample_name_list, len(sample_name_list), len(clone_dict_list)
	for (sample_name_X, clone_distribution_X) in zip(sample_name_list, clone_dict_list):
		for (sample_name_Y, clone_distribution_Y) in zip(sample_name_list, clone_dict_list):
			all_clone = set(clone_distribution_X.keys()) | set(clone_distribution_Y.keys())
			X_line, Y_line = [], []
			for clone_name in all_clone:
				try:
					X_line.append(math.log(float(clone_distribution_X[clone_name])))
				except KeyError:
					X_line.append(0)
				try:
					Y_line.append(math.log(float(clone_distribution_Y[clone_name])))
				except KeyError:
					Y_line.append(0)
			#if sample_name_X == sample_name_Y:
			#print zip(X_line, Y_line), len(zip(X_line, Y_line))
			#sys.exit(0)
			r = pearsonr(X_line, Y_line)
			print "Drawing %s and %s"%(sample_name_X, sample_name_Y)
			scatter_hist2(prj_folder, sample_name_X, sample_name_Y, X_line, Y_line)
			
			pickle_file = "%s/%s_%s_scatterinfo.txt"%(prj_folder, sample_name_X, sample_name_Y)
			pickle_file_handle = open(pickle_file, 'wb')
			dump_tuple = (X_line, Y_line, r)
			pickle.dump(dump_tuple, pickle_file_handle)
			pickle_file_handle.close()
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
