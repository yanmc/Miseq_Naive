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


def batch_iterator(iterator, batch_size) :
    """Returns lists of length batch_size.

    This can be used on any iterator, for example to batch up
    SeqRecord objects from Bio.SeqIO.parse(...), or to batch
    Alignment objects from Bio.AlignIO.parse(...), or simply
    lines from a file handle.

    This is a generator function, and it returns lists of the
    entries from the supplied iterator.  Each list will have
    batch_size entries, although the final list may be shorter.
    """
    entry = True #Make sure we loop once
    while entry :
        batch = []
        while len(batch) < batch_size :
            try :
                print help(iterator)
                entry = iterator.next()
            except StopIteration :
                entry = None
            if entry is None :
                #End of file
                break
            batch.append(entry)
        if batch :
            yield batch
def func(p,x):
	k,b=p
	return k*x+b
def error(p,x,y,s):
	#print s
	return func(p,x)-y

def main_v1():
	help(pd.read_csv)
	#sys.exit(0)
	mutation_patterns_files = glob.glob('%s/%s_*_mutation_patterns.txt'%(prj_tree.data, prj_name))
	mean_list = []
	for mutation_patterns_file in mutation_patterns_files:
		print mutation_patterns_file
		ref_seq_id_name = mutation_patterns_file.split("_")[-3]
		print ref_seq_id_name
		#if ref_seq_id_name == "IGHV1-18":
		mutation_patterns_file = open(mutation_patterns_file, "rU") 
		if os.fstat(mutation_patterns_file.fileno()).st_size:

			mutation_patterns_reader = np.loadtxt(mutation_patterns_file)
			if mutation_patterns_reader.ndim == 2:
				#print mutation_patterns_reader
				#print len(mutation_patterns_reader)
				#print mutation_patterns_reader.shape, mutation_patterns_reader.ndim
				line_numbers = np.sum(mutation_patterns_reader, axis=1)
				#print line_numbers
				line_number_percent_list = []
			elif mutation_patterns_reader.ndim == 1:
				#print mutation_patterns_reader
				mutation_patterns_reader = mutation_patterns_reader.reshape(1,len(mutation_patterns_reader))
				#print mutation_patterns_reader
				#print len(mutation_patterns_reader)
				#print mutation_patterns_reader.shape, mutation_patterns_reader.ndim
				line_numbers = np.sum(mutation_patterns_reader, axis=1)
				#print line_numbers
			for line_index, line in enumerate(mutation_patterns_reader):
				for row_index, row in enumerate(line):
					if line_numbers[line_index] != 0:
						mutation_patterns_reader[line_index, row_index] = mutation_patterns_reader[line_index, row_index] / line_numbers[line_index] * 100
			#print mutation_patterns_reader, mutation_patterns_reader.shape
			mutation_patterns_reader = copy.deepcopy(mutation_patterns_reader[ : 10])
			mutation_patterns_reader = mutation_patterns_reader.T
			#print mutation_patterns_reader, mutation_patterns_reader.shape
			mutation_patterns_reader_mean = np.mean(mutation_patterns_reader)
			mean_list.append(mutation_patterns_reader_mean)
			if mutation_patterns_reader.shape[1] == 1:
				plt.figure(figsize=(8,6))
				for line in mutation_patterns_reader:
					print Xi, Yi, type(Xi), type(Yi), len(Xi), len(Yi)
					Yi = line
					Xi = np.arange(1,len(line)+1)

					plt.scatter(Xi,Yi,color="red",linewidth=3)

				plt.legend()
				#plt.show()
			else:
				plt.figure(figsize=(8,6))
				for line in mutation_patterns_reader:
					Yi = line
					Xi = np.arange(1,len(line)+1)
					print Xi, Yi, type(Xi), type(Yi), len(Xi), len(Yi)
					p0=[10,2]
					s="Test the number of iteration"
					Para=leastsq(error,p0,args=(Xi,Yi,s))
					k,b=Para[0]

					if b >= 12.5 :
						print"k=",k,'\n',"b=",b

						plt.scatter(Xi,Yi,color="red",label="Sample Point",linewidth=3)
						plt.plot(Xi,Yi,color="red",linewidth=2)
						x=np.linspace(0,10,1000)
						y=k*x+b
						plt.plot(x,y,color="orange",label="Fitting Line",linewidth=2)
					else:
						plt.scatter(Xi,Yi,color="blue",linewidth=3)
						plt.plot(Xi,Yi,color="blue",linewidth=2)
				plt.legend()
				plt.show()
				#sys.exit(0)


			'''
			print mutation_patterns_reader
			index, sum_index = 0, 0
			while index <= len(mutation_patterns_reader):

				index += 10
				batch = mutation_patterns_reader[sum_index : index]
				sum_index = index
				print batch
			'''	
		else:
			print "yes"
	#mean_mutation = float(sum(mean_list))/float(len(mean_list))
	#print mean_mutation
		
		#for line_number in line_numbers:
		#	line_number_percent = float(line_number)/float(total_number) * 100
		#	line_number_percent_list.append(line_number_percent)
		#print len(line_number_percent_list), len(columns)
		
		

def main():
	help(pd.read_csv)
	#sys.exit(0)
	mutation_patterns_files = glob.glob('%s/%s_*_mutation_patterns.txt'%(prj_tree.data, prj_name))
	mean_list = []
	for mutation_patterns_file in mutation_patterns_files:
		print mutation_patterns_file
		ref_seq_id_name = mutation_patterns_file.split("_")[-3]
		print ref_seq_id_name
		#if ref_seq_id_name == "IGHV1-18":
		mutation_patterns_file = open(mutation_patterns_file, "rU") 
		if os.fstat(mutation_patterns_file.fileno()).st_size:

			mutation_patterns_reader = np.loadtxt(mutation_patterns_file)
			if mutation_patterns_reader.ndim == 2:
				#print mutation_patterns_reader
				#print len(mutation_patterns_reader)
				#print mutation_patterns_reader.shape, mutation_patterns_reader.ndim
				line_numbers = np.sum(mutation_patterns_reader, axis=1)
				#print line_numbers
				line_number_percent_list = []
			elif mutation_patterns_reader.ndim == 1:
				#print mutation_patterns_reader
				mutation_patterns_reader = mutation_patterns_reader.reshape(1,len(mutation_patterns_reader))
				#print mutation_patterns_reader
				#print len(mutation_patterns_reader)
				#print mutation_patterns_reader.shape, mutation_patterns_reader.ndim
				line_numbers = np.sum(mutation_patterns_reader, axis=1)
				#print line_numbers
			for line_index, line in enumerate(mutation_patterns_reader):
				for row_index, row in enumerate(line):
					if line_numbers[line_index] != 0:
						mutation_patterns_reader[line_index, row_index] = mutation_patterns_reader[line_index, row_index] / line_numbers[line_index] * 100
			#print mutation_patterns_reader, mutation_patterns_reader.shape
			mutation_patterns_reader = copy.deepcopy(mutation_patterns_reader[ : 10])
			mutation_patterns_reader = mutation_patterns_reader.T

			#print mutation_patterns_reader, mutation_patterns_reader.shape
			mutation_patterns_reader_mean = np.mean(mutation_patterns_reader)
			mean_list.append(mutation_patterns_reader_mean)
			if mutation_patterns_reader.shape[1] == 1:
				plt.figure(figsize=(8,6))
				for line in mutation_patterns_reader:
					print Xi, Yi, type(Xi), type(Yi), len(Xi), len(Yi)
					Yi = line
					Xi = np.arange(1,len(line)+1)

					plt.scatter(Xi,Yi,color="red",linewidth=3)

				plt.legend()
				#plt.show()
			else:
				left, width = .45, .5
				bottom, height = .5, .5
				right = left + width
				top = bottom + height
				fig = plt.figure(figsize=(8,6))
				ax = fig.add_subplot(111)
				line_flag = 1
				for index, line in enumerate(mutation_patterns_reader):
					Yi = line
					Xi = np.arange(1,len(line)+1)
					#print Xi, Yi, type(Xi), type(Yi), len(Xi), len(Yi)
					X=sm.add_constant(Xi) 
					est=sm.OLS(Yi,X)
					
					est=est.fit()
					if est.params[0] >= 12.5:
						
						print est.summary()
						print est.params
						print est.tvalues
						print est.pvalues
						#print help(sm.regression.linear_model.OLSResults)
						#sys,exit(0)
						X_prime=np.linspace(X.min(), X.max(),100)[:,np.newaxis]
						X_prime=sm.add_constant(X_prime)
						y_hat=est.predict(X_prime)
						plt.plot(Xi, Yi, 'g', alpha=0.9, linewidth=2)
						plt.scatter(Xi, Yi, c = 'g', marker = "o", alpha=0.9, linewidth=2)
						
						plt.plot(X_prime[:,1], y_hat, 'r')
						ax.text(right, top - (0.05 * line_flag), "Nt postition: " + str(index+1), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
						ax.text(right, top - (0.05 * line_flag)-0.05, "P-value: " + str(round(est.pvalues[0], 3)), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
						line_flag += 2
						
					else:
						plt.plot(Xi, Yi, 'b', alpha=0.3, linewidth=2)
						plt.scatter(Xi, Yi, c = 'b', marker = "o", alpha=0.3, linewidth=2)
						
				plt.xlabel("Mutation Count (Sequence)")
				plt.ylabel("Mutation Freq. (position)")
				plt.ylim(ymin=0)
				#plt.ylim(0, 100)
				plt.xlim(1, 10)
				plt.title(ref_seq_id_name)
				fig.savefig("%s/%s_%s_mutation_patterns.backup.png"%(prj_tree.figure, prj_name, ref_seq_id_name), dpi = 300)
				plt.close(fig)
				#sys.exit(0)


			'''
			print mutation_patterns_reader
			index, sum_index = 0, 0
			while index <= len(mutation_patterns_reader):

				index += 10
				batch = mutation_patterns_reader[sum_index : index]
				sum_index = index
				print batch
			'''	
		else:
			print "yes"
	#mean_mutation = float(sum(mean_list))/float(len(mean_list))
	#print mean_mutation
		
		#for line_number in line_numbers:
		#	line_number_percent = float(line_number)/float(total_number) * 100
		#	line_number_percent_list.append(line_number_percent)
		#print len(line_number_percent_list), len(columns)

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