#!/usr/bin/env python
# encoding: utf-8
"""
misc_detect_primer_efficiency.py

Created by Mingchen on 2015-01-19.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.

Usage: misc-get-pattern-loc.py -i infile -r referencefile -o outfile


"""
import sys, os, csv, re, glob, copy, subprocess, time, multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Process, Manager
from common_info import *
from mytools import *
from misc_prepare_pbs import *
from collections import Counter 
try:
    import cPickle as pickle
except ImportError:
    import pickle


def count_primer(prj_tree, prj_name,mapdict):
	infile = glob.glob("%s/*.fasta" %(prj_tree.origin))[0]
	writer = csv.writer(open("%s/%s_primer_distribution.txt"%(prj_tree.data, prj_name),'w'),delimiter ='\t')
	writer2 = csv.writer(open("%s/%s_primer_in_seq.txt"%(prj_tree.data, prj_name),'w'),delimiter ='\t')
	writer3 = csv.writer(open("%s/%s_primer_for_plot.txt"%(prj_tree.data, prj_name),'w'),delimiter ='\t')
	writer.writerow(['primer','primer_seq','seq_num','seq_id','gene_num','gene_id'])
	writer2.writerow(['seq_id','primer_num','primer_id'])

	recordlist, noprimer, pprimer, rprimer = [], 0, 0, 0
	primerdict = {}
	
	for key in MULTIPLEX_PRIMER.keys():
		primerdict.setdefault(key,[])
	
	for record in SeqIO.parse(infile,'fasta'):
		recordlist.append(record.id)
		myseq = record.seq
		seqlength = len(str(myseq))
		flag = []
		for primers in MULTIPLEX_PRIMER:
			starts = str(myseq).find(MULTIPLEX_PRIMER[primers])
			if starts != -1 :
				primerdict[primers].append(record.id)
				flag.append(primers)
			else:
				starts = str(myseq.reverse_complement()).find(MULTIPLEX_PRIMER[primers])
				if starts != -1:
					primerdict[primers].append(record.id)
					flag.append(primers)

		if len(flag) >= 1:				
			writer2.writerow([record.id,len(flag)] + flag) 
		if len(flag) == 0:
			noprimer += 1

	print '%d reads have no primer !'%noprimer
	titleflag = 0
	for key in primerdict.keys():
		genes = []
		for iterms in primerdict[key]:
			#print iterms
			genes.append(mapdict[iterms])
		mycount = Counter(genes)
		writer.writerow([key,MULTIPLEX_PRIMER[key],len(primerdict[key]),primerdict[key],len(genes),genes])
		
		plot_dict = {}
		for family in HUMAN_GERMLINE_V_ALLELS.keys():
			if 'IGHV' in family:
				for allels in HUMAN_GERMLINE_V_ALLELS[family]:
					if allels in mycount.keys():
						plot_dict[allels] = mycount[allels]
						print mycount[allels]
					else:
						plot_dict[allels] = 0
		mylist,writelist = [],[]
		for mykey in sorted(plot_dict.keys()):
			mylist.append(plot_dict[mykey])
			#keylist.append(mykey)
		if titleflag == 0:
			writer3.writerow(sorted(plot_dict.keys()))
			titleflag += 1
		else:
			pass
		mylist.insert(0,sum(mylist))
		mylist.insert(0,key)	
		writer3.writerow(mylist)
	
	
def get_align(prj_tree, prj_name):
	align_info = glob.glob("%s/%s_get_recombanation_info.txt" %(prj_tree.igblast_data,prj_name))[0] 	
	aligns = csv.reader(open(align_info,'r'),delimiter = '\t')
	mapdict = {}
	for line in aligns:
		if len(line) > 1:
			mapdict[line[0].strip()] = line[1].split(',')[0]
		elif len(line) == 1:
			mapdict[line[0].strip()] = ''
	return mapdict	

def main():
	mapdict = get_align(prj_tree, prj_name)
	count_primer(prj_tree, prj_name,mapdict)
if __name__=='__main__':
	prj_folder = os.getcwd()
	prj_tree = ProjectFolders(prj_folder)
	prj_name = fullpath2last_folder(prj_tree.home)
	
	main()
