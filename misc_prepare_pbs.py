#!/usr/bin/env python
# encoding: utf-8
"""
misc_prepare_pbs.py

Created by Mingchen on 2015-01-19.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.

Usage: misc-get-pattern-loc.py -i infile -r referencefile -o outfile


"""
import sys, os, csv, re, glob, copy, subprocess, time, multiprocessing
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, Process, Manager
#from bsub import bsub
from mytools import *
from misc_prepare_pbs import *
try:
    import cPickle as pickle
except ImportError:
    import pickle

def processing_jobs(job):
	cmd = '%s  %s'%('bash', job)
	head, tail 	= os.path.splitext(job)	
	p=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
def bsub_jobs(job):
	cmd = '%s <  %s'%('bsub', job)
	head, tail 	= os.path.splitext(job)	
	p=subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)  
	while True:
		buff = p.stdout.readline()
		if ">" in buff and "<" in buff:
	 		igblast_job_id = buff.split(">")[0].split("<")[-1]
		if buff == '' and p.poll() != None:
			break
	return igblast_job_id

def check_jobs_done(prj_name, prj_tree, app, igblast_job_ids):
	log_file = "%s/%s_pbs.log"%(prj_tree.logs, app)
	log_file_handle = open(log_file, "w")
	for index, igblast_job_id in enumerate(igblast_job_ids):
		while os.path.exists("%s/output_%s"%(prj_tree.jobs, igblast_job_id)) == False:
			log_file_handle.write("Waiting for job_%s...%s\n"%(igblast_job_id,time.ctime()))
			time.sleep(10)
		
		errput, output = "%s/errput_%s"%(prj_tree.jobs, igblast_job_id), "%s/output_%s"%(prj_tree.jobs, igblast_job_id)
		IgBLAST_log = True
		while IgBLAST_log:
			output_log = open(output, "rU")
			for line in output_log.readlines():
				if line.replace('\n','') == "Successfully completed.":
					log_file_handle.write("job_%s Successfully completed.\n"%igblast_job_id)
					IgBLAST_log = False
					#os.system("mv %s %s"%(errput, prj_tree.logs))
					#os.system("mv %s %s"%(output, prj_tree.logs))
					#os.system("rm %s"%(errput))
					os.system("rm %s"%(output))
			output_log.close()

		#igblast_job_ids.pop(index)

def prepare_cdhit_nucle_pbs(prj_name, prj_tree, fasta_file, round_index):
	head, tail 	= os.path.splitext(fasta_file)
	barcode 	= head.split("/")[-1].split("_")[4]
	handle = open("%s/prepare_cdhit_nucle_pbs_%s_%s.sh" %(prj_tree.jobs, barcode, round_index), "w")
	handle.write("#!/bin/bash\n")
	handle.write("#BSUB -J %s_%s_%s\n" %("prepare_cdhit_nucle_pbs", barcode, round_index))
	handle.write("#BSUB -n 1\n")
	#handle.write("#BSUB -n %s\n"%(infile_number*4))
	handle.write("#BSUB -R %s\n"%("\"span[ptile=1]\""))
	handle.write("#BSUB -o %s/output_%%%s\n"%(prj_tree.jobs, "J"))
	handle.write("#BSUB -e %s/errput_%%%s\n"%(prj_tree.jobs, "J"))
	handle.write("#BSUB -q cpu\n")
	handle.write("cd-hit-est -i %s -o %s_%s -c 0.95 -n 10 -d 0 -M 0 -T 8 -g 1 -p 1 &"%(fasta_file, head, round_index))
	handle.close()
def prepare_justfy_primer_and_group_pbs(prj_name, prj_tree, pickle_file):
	head, tail 	= os.path.splitext(pickle_file)
	barcode 	= head.split("/")[-1].split("_")[-1]
	handle = open("%s/justfy_primer_and_group_%s.sh" %(prj_tree.jobs, barcode), "w")
	handle.write("#!/bin/bash\n")
	handle.write("#BSUB -J %s_%s\n" %("justfy_primer_and_group", barcode))
	handle.write("#BSUB -n 1\n")
	#handle.write("#BSUB -n %s\n"%(infile_number*4))
	handle.write("#BSUB -R %s\n"%("\"span[ptile=1]\""))
	handle.write("#BSUB -o %s/output_%%%s\n"%(prj_tree.jobs, "J"))
	handle.write("#BSUB -e %s/errput_%%%s\n"%(prj_tree.jobs, "J"))
	handle.write("#BSUB -q cpu\n")
	handle.write("python ./2.1-justfy-primer-and-group.py -i %s &"%(pickle_file))
	handle.close()
def prepare_clustal_jobs_normal(prj_name, prj_tree, UMI_length, group_type):
	clustal_fastas = glob.glob("%s/%s_*_cut_berfore_UMI%s_%s.fasta"%(prj_tree.clustal_fasta, prj_name, UMI_length, group_type))
	for infile in clustal_fastas:
		head, tail 	= os.path.splitext(infile)
		fname		= head.split("/")[-1]
		handle = open("%s/clustal_%s.sh" %(prj_tree.jobs, fname), "w")
		handle.write("#!/bin/bash\n")
		handle.write("#BSUB -J %s_%s\n" %(prj_name, fname))
		handle.write("#BSUB -n 1\n")
		#handle.write("#BSUB -n %s\n"%(infile_number*4))
		handle.write("#BSUB -R %s\n"%("\"span[ptile=1]\""))
		handle.write("#BSUB -o %s/output_%%%s\n"%(prj_tree.jobs, "J"))
		handle.write("#BSUB -e %s/errput_%%%s\n"%(prj_tree.jobs, "J"))
		handle.write("#BSUB -q cpu\n")
		handle.write("/zzh_gpfs/apps/clustalw-2.1-linux-x86_64-libcppstatic/clustalw2 -infile=%s  -ITERATION=ALIGNMENT&"%(infile))
		handle.close()

def prepare_IgBLAST_jobs(prj_name, prj_tree):
	"""
	prepare files for each of the fasta files todo IgBLAST
	"""
	
	# load all fasta files in split folder
	infiles = glob.glob("%s/%s_*.fasta" %(prj_tree.split, prj_name))
	infile_number = len(infiles)
	
	for infile in infiles:
		head, tail 	= os.path.splitext(infile)
		f_ind 		= head.split("_")[ -1 ]
		handle = open("%s/IgBLAST_%s.sh" %(prj_tree.jobs,f_ind), "w")
		handle.write("#!/bin/bash\n")
		handle.write("#BSUB -J %s_%s\n" %(prj_name,f_ind))
		handle.write("#BSUB -n 1\n")
		#handle.write("#BSUB -n %s\n"%(infile_number*4))
		handle.write("#BSUB -R %s\n"%("\"span[ptile=1]\""))
		handle.write("#BSUB -o %s/output_%%%s\n"%(prj_tree.jobs, "J"))
		handle.write("#BSUB -e %s/errput_%%%s\n"%(prj_tree.jobs, "J"))
		handle.write("#BSUB -q cpu\n")
		handle.write("igblastn -germline_db_V ./Igblast_database/20150429-human-gl-v -germline_db_J \
		./Igblast_database/20150429-human-gl-j -germline_db_D ./Igblast_database/20150429-human-gl-d \
		-organism human -domain_system imgt -query %s -auxiliary_data optional_file/human_gl.aux \
		-outfmt '7 qseqid sseqid pident length mismatch gapopen gaps qstart qend sstart send evalue \
		bitscore qlen slen qseq sseq score frames qframe sframe positive ppos btop staxids stitle \
		sstrand qcovs qcovhsp' -num_alignments_V 3 -num_alignments_D 3 -num_alignments_J 3 -out \
		%s/IgBLAST_result_%s.txt &"%(infile, prj_tree.igblast_data, f_ind))
		handle.close()

def main():
	print "This is a module!"
if __name__=='__main__':
	main()