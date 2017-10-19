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

def main():
	#sys.exit(0)
	
	germline_fasta = SeqIO.index("/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/Naive_IgM/Igblast_database/20150429-human-gl-vdj.fasta", "fasta")
	
	
	
	mutation_patterns_files = glob.glob('%s/%s_*_mutation_patterns.txt'%(prj_tree.data, prj_name))
	mean_list = []
	ploy_gene_position_dict = {}
	for mutation_patterns_file in mutation_patterns_files:
		print mutation_patterns_file
		ref_seq_id_name = mutation_patterns_file.split("_")[-3]
		print ref_seq_id_name
		#if ref_seq_id_name == "IGHV1-18":
		mutation_patterns_file = open(mutation_patterns_file, "rU") 
		if os.fstat(mutation_patterns_file.fileno()).st_size:
			
			mutation_patterns_reader = np.loadtxt(mutation_patterns_file)
			print mutation_patterns_reader, len(mutation_patterns_reader)
			mutation_patterns_reader = copy.deepcopy(mutation_patterns_reader[ : 10])
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
			
			line_numbers = [x/(index+1) for (index, x) in enumerate(line_numbers)]
			third_quartile = scist.scoreatpercentile(line_numbers,75)
			interquartile_range = scist.scoreatpercentile(line_numbers,75) - scist.scoreatpercentile(line_numbers,25)
			#print mutation_patterns_reader, len(mutation_patterns_reader)
			print third_quartile, interquartile_range, line_numbers
			start_bin = 0
			for bin_index in sorted(range(1,6), reverse = True):
				#print bin_index, mutation_patterns_reader[bin_index]
				try: 
					bin_numbers =  line_numbers[bin_index]
					print bin_numbers, third_quartile, interquartile_range, ( bin_numbers - third_quartile) - (1.5 * interquartile_range)
					if ( bin_numbers - third_quartile) > (1.5 * interquartile_range):#0.4 Naive :0.3
						print "yes"
						print bin_numbers, third_quartile, interquartile_range, ( bin_numbers - third_quartile) - (1.5 * interquartile_range)
						mutation_patterns_reader = copy.deepcopy(mutation_patterns_reader[bin_index:])
						start_bin = bin_index
						break
				except IndexError:
					continue
			#if ref_seq_id_name == "IGHV1-69":
			#	sys.exit(0)
			for line_index, line in enumerate(mutation_patterns_reader):
				for row_index, row in enumerate(line):
					if line_numbers[line_index] != 0:
						#print line_index, line_index + start_bin, line_numbers[line_index + start_bin], line_numbers
						mutation_patterns_reader[line_index, row_index] = mutation_patterns_reader[line_index, row_index] / line_numbers[line_index + start_bin] * 100
			#print mutation_patterns_reader, mutation_patterns_reader.shape
			
			mutation_patterns_reader = mutation_patterns_reader.T
			#print mutation_patterns_reader, mutation_patterns_reader.shape
			if mutation_patterns_reader.shape[1] == 1:
				ploy_gene_position_dict[ref_seq_id_name] = [0]
			else:
				ploy_gene_position_dict[ref_seq_id_name] = [0]
				left, width = .45, .5
				bottom, height = .5, .5
				right = left + width
				top = bottom + height
				fig = plt.figure(figsize=(8,6))
				ax = fig.add_subplot(111)
				line_flag = 1
				for index, line in enumerate(mutation_patterns_reader):
					Yi = line
					Xi = np.arange(start_bin + 1 , start_bin + len(line) + 1)
					#print Xi, Yi, type(Xi), type(Yi), len(Xi), len(Yi)
					X=sm.add_constant(Xi) 
					est=sm.OLS(Yi,X)
					est=est.fit()
					if est.params[0] >= 12.5:
						if est.pvalues[0] <= 0.05:
							ploy_gene_position_dict.setdefault(ref_seq_id_name, []).append(index+1)
							#print est.summary()
							print est.params[0]
							#print est.tvalues
							print est.pvalues[0]
							#print help(sm.regression.linear_model.OLSResults)
							#sys,exit(0)
							X_prime=np.linspace(X.min(), X.max(),100)[:,np.newaxis]
							X_prime=sm.add_constant(X_prime)
							y_hat=est.predict(X_prime)
							plt.plot(Xi, Yi, 'g', alpha=0.9, linewidth=2)
							plt.scatter(Xi, Yi, c = 'g', marker = "o", alpha=0.9, linewidth=2)

							plt.plot(X_prime[:,1], y_hat, 'r')
							ax.text(.2, top - (0.05 * line_flag), "Nt postition: " + str(index+1), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
							ax.text(.2, top - (0.05 * line_flag)-0.05, "P-value: " + str(round(est.pvalues[0], 3)), horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
							line_flag += 2
						
					else:
						plt.plot(Xi, Yi, 'b', alpha=0.3, linewidth=2)
						plt.scatter(Xi, Yi, c = 'b', marker = "o", alpha=0.3, linewidth=2)
						
				plt.xlabel("Mutation Count (Sequence)")
				plt.ylabel("Mutation Freq. (position)")
				plt.ylim(ymin=0)
				plt.ylim(0, 100)
				#plt.xlim(start_bin , start_bin + len(line))
				plt.xlim(0 , 10)
				#plt.xticks(range(0,len(line),1), ('0','1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))
				plt.title(ref_seq_id_name)
				fig.savefig("%s/%s_%s_mutation_patterns0.4.png"%(prj_tree.figure, prj_name, ref_seq_id_name), dpi = 300)
				plt.close(fig)
		else:
			print "yes"
			ploy_gene_position_dict[ref_seq_id_name] = [0]
			
	#output 
	print ploy_gene_position_dict
	pickle_file = '%s/%s_ploy_gene_position_dict_dump'%(prj_tree.tmp, prj_name)
	pickle_file_handle = open(pickle_file, 'wb')
	dump_tuple = (ploy_gene_position_dict)
	pickle.dump(dump_tuple, pickle_file_handle)
	pickle_file_handle.close()
	#Step 2: get PLOY position nucl and plot stack bar
	for germline_type in ('V','J'):
		pickle_file = '%s/%s_gene_usage_info_dump_%s_%s%s'%(prj_tree.tmp, prj_name, germline_type, chain_type, pic_type)
		f_IgH = open(pickle_file, 'rb')
		pickle_tuple_IgH = pickle.load(f_IgH)
		columns_IgH, geneusage_data_IgH, gene_usage_ids_IgH = pickle_tuple_IgH[0], pickle_tuple_IgH[1], pickle_tuple_IgH[2]
		f_IgH.close()
		geneusage_dict = {}
		for i in range(len( pickle_tuple_IgH[0])):
			geneusage_dict.setdefault(columns_IgH[i].split("*")[0], []).append(sum(geneusage_data_IgH[i]))
		for (key, value) in geneusage_dict.items():
			geneusage_dict[key] = sum(value)
		result_file_name = open('%s/%s_%s_%s_right_allele_usage%s.txt'%(prj_tree.data, prj_name, chain_type, germline_type, pic_type),"rU")
		result_file = csv.reader(result_file_name, delimiter="\t")
		max_freq_allele_dict = {}
		for line in result_file:
			max_freq_allele_dict[line[0]] = line[0] + "*" + line[1]
		print max_freq_allele_dict
		
		pickle_files = glob.glob('%s/%s_mutation_pattrens_dump_%s_%s_*'%(prj_tree.tmp, prj_name, germline_type,chain_type))
		for pickle_file in pickle_files:
			
			print pickle_file
			ref_seq_id_name = pickle_file.split('_')[-1]
			outfile = open('%s/%s_mutation_spectrum_%s_%s_%s_all'%(prj_tree.data, prj_name, germline_type, chain_type, ref_seq_id_name), "w")
			writer= csv.writer(outfile, delimiter= "\t")
			pickle_file_handle = open(pickle_file, 'rb')
			pickle_tuple = pickle.load(pickle_file_handle)
			mutation_patterns_group = pickle_tuple
			#print mutation_patterns_group
			pickle_file_handle.close()
			try:
				max_freq_allele = max_freq_allele_dict[ref_seq_id_name]
			except KeyError:
				continue
			germline_fasta_seq = germline_fasta[max_freq_allele].seq.upper()
			mutation_spectrum_array = np.zeros((len(germline_fasta_seq), 5))
			for (key, value) in mutation_patterns_group.items():
				#if key <= 10:
				for record in value:
					read_id = record[0]
					position_records = record[1]
					for position_record_item in position_records:
						position = position_record_item[0]
						ref_nucl = position_record_item[1]
						query_nucl = position_record_item[2]
						if query_nucl == "A":
							mutation_spectrum_array[position-1][0] += 1
						elif query_nucl == "C":
							mutation_spectrum_array[position-1][1] += 1
						elif query_nucl == "-":
							mutation_spectrum_array[position-1][2] += 1
						elif query_nucl == "T":
							mutation_spectrum_array[position-1][3] += 1
						elif query_nucl == "G":
							mutation_spectrum_array[position-1][4] += 1
			for index, nucl in enumerate(germline_fasta_seq):
				if nucl == "A":
					print "mutation_spectrum_array[index][0]", mutation_spectrum_array[index][0]
					#mutation_spectrum_array[index][0] =  0#-sum(mutation_spectrum_array[index])
				if nucl == "C":
					print "mutation_spectrum_array[index][1]", mutation_spectrum_array[index][1]
					#mutation_spectrum_array[index][1] =  0#-sum(mutation_spectrum_array[index])
				if nucl == "T":
					print "mutation_spectrum_array[index][3]", mutation_spectrum_array[index][3]
					#mutation_spectrum_array[index][3] =  0#-sum(mutation_spectrum_array[index])
				if nucl == "G":
					print "mutation_spectrum_array[index][4]", mutation_spectrum_array[index][4]
					#mutation_spectrum_array[index][4] =  0#-sum(mutation_spectrum_array[index])
			for index, line in enumerate(mutation_spectrum_array):			
				writer.writerow([germline_fasta_seq[index]] + list(line))
			outfile.close()
			#Plot mutation spectrum
		 	d_colors = ['#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777']
			#gap = 0.2
			SBG = StackedBarGrapher()
			fig = plt.figure()
			ax5 = fig.add_subplot(111)
			SBG.stackedBarPlot(ax5,mutation_spectrum_array,d_colors,edgeCols=d_colors,ylabel = 'Number of reads')
			plt.title( "%s %s"%(ref_seq_id_name, geneusage_dict[ref_seq_id_name]))
			#plt.xticks(range(0,10,1), ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))
			
			plt.savefig('%s/%s_%s_%s_%s_nutation_spectrum_all.png'%(prj_tree.figure, prj_name, chain_type, germline_type, ref_seq_id_name ), dpi=300)
			del fig
			plt.close()
			
			
			# ploy_nucl_position_draw
			poly_nucl_positions = ploy_gene_position_dict[ref_seq_id_name]
			if len(poly_nucl_positions) > 1:
				for poly_nucl_position in poly_nucl_positions[1:]:
					poly_nucl_position_array = np.zeros((10, 4))
					
					for (key, value) in mutation_patterns_group.items():
						if key <= 10:
							for record in value:
								read_id = record[0]
								position_records = record[1]
								for position_record_item in position_records:
									position = position_record_item[0]
									if position == poly_nucl_position:
										ref_nucl = position_record_item[1]
										query_nucl = position_record_item[2]
										if query_nucl == "A":
											poly_nucl_position_array[key-1][0] += 1
										elif query_nucl == "T":
											poly_nucl_position_array[key-1][1] += 1
										elif query_nucl == "C":
											poly_nucl_position_array[key-1][2] += 1
										elif query_nucl == "G":
											poly_nucl_position_array[key-1][3] += 1

					#print poly_nucl_position_array
					#d_colors = ['#2166ac', '#fee090', '#fdbb84', '#fc8d59', '#e34a33', '#b30000', '#777777']
					d_colors = ['red', 'yellow', 'blue', 'green']
					gap = 0.2
					SBG = StackedBarGrapher()
					fig = plt.figure()
					ax5 = fig.add_subplot(111)
					SBG.stackedBarPlot(ax5,poly_nucl_position_array,d_colors,edgeCols=['#000000']*7,ylabel = 'Number of reads',gap=gap)
					plt.title( "%s %s %s"%(ref_seq_id_name, poly_nucl_position, ref_nucl))
					plt.xticks(range(0,10,1), ('1', '2', '3', '4', '5', '6', '7', '8', '9', '10'))
					plt.savefig('%s/%s_%s_%s_%s_position%s_all.png'%(prj_tree.figure, prj_name, chain_type, germline_type, ref_seq_id_name, poly_nucl_position ), dpi=300)
					del fig
					plt.close()
		
if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	start = time.time()
	pic_type = ""#"_unique"
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