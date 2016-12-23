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
	certain_germline = ["IGHV1-58*01", "IGHV1-69*01", "IGHV3-30*04", "IGHV4-28*03", "IGHV3-7*01", "IGHV4-59*08", "IGHV6-1*01", "IGHV7-4-1*01", "IGHV6-1*01", "IGHV4-61*01", "IGHV3-48*01", "IGHV1-2*01"]
	germline_fasta = SeqIO.index("/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/Naive_IgM/Igblast_database/20150429-human-gl-vdj.fasta", "fasta")
	
	pickle_file = '%s/%s_ploy_gene_position_dict_dump'%(prj_tree.tmp, prj_name)
	ref_seq_id_name = pickle_file.split('_')[-1]
	pickle_file_handle = open(pickle_file, 'rb')
	pickle_tuple = pickle.load(pickle_file_handle)
	ploy_gene_position_dict = pickle_tuple
	#print ploy_gene_position_dict
	ploy_gene_position_dict["IGHV3-30"] = [0,98,99,180,293]
	ploy_gene_position_dict["IGHV3-7"] = [0,276]
	ploy_gene_position_dict["IGHV4-28"] = [0,173,242,295,296]
	#print ploy_gene_position_dict
	
	for germline_type in ('V','J'):	
		
		pickle_files = glob.glob('%s/%s_mutation_pattrens_dump_%s_%s_*'%(prj_tree.tmp, prj_name, germline_type,chain_type))
		for pickle_file in pickle_files:

			print pickle_file
			ref_seq_id_name = pickle_file.split('_')[-1]
			pickle_file_handle = open(pickle_file, 'rb')
			pickle_tuple = pickle.load(pickle_file_handle)
			mutation_patterns_group = pickle_tuple
			
			poly_nucl_positions = ploy_gene_position_dict[ref_seq_id_name]
			
			reads_poly_comb = []
			for (key, value) in mutation_patterns_group.items():
				#print key
				if key <= 10:
					for record in value:
						read_poly_comb = ''
						read_id = record[0]
						
						position_records = record[1]
						for position_record_item in position_records:
							position = position_record_item[0]
							if len(poly_nucl_positions) > 1:
								for poly_nucl_position in poly_nucl_positions[1:]:
									if position == poly_nucl_position:
										ref_nucl = position_record_item[1]
										query_nucl = position_record_item[2]
										read_poly_comb = read_poly_comb + str(position) +'_'+ query_nucl + '_'
						if read_poly_comb != '':
							reads_poly_comb.append(read_poly_comb)
			#print ref_seq_id_name, list_counter(reads_poly_comb)
			for allele in certain_germline:
				gene_name = allele.split("*")[0]
				if ref_seq_id_name == gene_name:
					certain_seq = str(germline_fasta[allele].seq)
					for (key, value) in  list_counter(reads_poly_comb).most_common(3):
						novel_seq = certain_seq
						print  key, value
						print "It's a novel gene"
						nucls = key.split("_")[1::2]
						poss = key.split("_")[:-1:2]
						#print poss, nucls
						for (x, y) in zip(poss, nucls):
							print "Position: %s, nucleotide: %s"%(x, y)
							novel_seq = novel_seq[ : int(x)-1] + y + certain_seq[int(x) : ]
							
						#print certain_seq
						#print novel_seq
						#print str(germline_fasta['IGHV1-58*02'].seq)
						
						for record in germline_fasta:
							if novel_seq.upper() == str(germline_fasta[record].seq).upper() or novel_seq.upper() in str(germline_fasta[record].seq).upper() or str(germline_fasta[record].seq).upper() in novel_seq.upper():
								#print novel_seq.upper()
								#print str(germline_fasta[record].seq).upper()
								#print ref_seq_id_name,allele, poss, nucls, record
								print "It's %s"%record
								
			
	
if __name__ == '__main__':

	pool_size = multiprocessing.cpu_count()
	# create 1st and 2nd subfolders
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	start = time.time()
	pic_type = ""
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