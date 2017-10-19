#!/usr/bin/env python
# encoding: utf-8
"""
misc-choose-spike-in.py

Created by Mingchen on 2015-05-15.
Copyright (c) 2015 __MyCompanyName__. All rights reserved.
"""

import sys
import os
import csv
import re
import glob
import copy
import subprocess
import time
import multiprocessing
print sys.path
import traceback
import tempfile
import numpy as np
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
try:
    import cPickle as pickle
except ImportError:
    import pickle


def get_all_gene_name(germline_gene_list):
    return sorted(set([x.split('*')[0] for x in germline_gene_list]))


def get_germline_gene(germline_type, chain_type):
    if germline_type == 'V' and chain_type == "H":
        germline_gene_list = HUMAN_GERMLINE['HUMANIGHV']
    elif germline_type == 'J' and chain_type == "H":
        germline_gene_list = HUMAN_GERMLINE['HUMANIGHJ']

    elif germline_type == 'V' and chain_type == "K":
        germline_gene_list = HUMAN_GERMLINE['HUMANIGKV']
    elif germline_type == 'J' and chain_type == "K":
        germline_gene_list = HUMAN_GERMLINE['HUMANIGKJ']

    elif germline_type == 'V' and chain_type == "L":
        germline_gene_list = HUMAN_GERMLINE['HUMANIGLV']
    elif germline_type == 'J' and chain_type == "L":
        germline_gene_list = HUMAN_GERMLINE['HUMANIGLJ']
    return germline_gene_list


def load_assignment_dict(IGBLAST_assignment_file):
    assignment_dict = {}
    infile = open(IGBLAST_assignment_file, "rU")
    reader = csv.reader(infile, delimiter="\t")
    for index, line in enumerate(reader):
        assign_result = MyAlignment(line)
        assignment_dict.setdefault(assign_result.qid,
                                   []).append([assign_result.assign_type,
                                               assign_result.sstart,
                                               assign_result.qseq,
                                               assign_result.sseq,
                                               assign_result.sid])
    infile.close()
    return assignment_dict


def main():
    #germline_fasta = SeqIO.index("/zzh_gpfs02/yanmingchen/HJT-PGM/Naive/Naive_IgM/Igblast_database/20150429-human-gl-vdj.fasta", "fasta")
    #unique_real_reads_fasta = SeqIO.index("%s/%s_%s_real_reads_Variable_region.fasta"%(prj_tree.reads, prj_name, chain_type), "fasta")
    #IGBLAST_assignment_file = "%s/%s_get_assignment_info.txt"%(prj_tree.igblast_data, prj_name)
    #assignment_dict = load_assignment_dict(IGBLAST_assignment_file)
    print sys.path
    primer_efficiency_file = open(
        "%s/%s_primer_distribution.txt" %
        (prj_tree.data, prj_name), "rU")
    primer_efficiency_handle = csv.reader(
        primer_efficiency_file, delimiter="\t")
    for line in primer_efficiency_handle:
        print line


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
    print prj_name, end - start
    print "Finished"
