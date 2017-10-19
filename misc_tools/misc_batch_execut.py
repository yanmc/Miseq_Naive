#!/usr/bin/env python
# encoding: utf-8
"""
3.0.py

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
import traceback
import tempfile
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
        # and( "K" in file_name or "L" in file_name):
        if os.path.isdir(file_name) :#and "TCR" in file_name:
            os.chdir(file_name)
            num = file_name[-3:]
            print os.getcwd()  # , os.system("which 2.0.py")
            
                
            #os.system("bsub -n 1 -q zzh -e err.10 -o out.10 -P %s \"misc_getR1R2_reads.py\""%("R1R2" ))
            #os.system("mv ./origin/output/M1_atleast-2.fasta ./origin/%s.assembled.fasta"%file_name)
            
            #os.system("bsub -n 2 -q zzh -e err.6 -o out.6 -P %s \"misc_del_mutibarcodeprimer.py\""%("del" ))
            #os.system("bsub -n 1 -R \"span[hosts=1]\" -q zzh -e err.1 -o out.1 \"1.0.py\"")
            #os.system("grep \"GCGTAGTACTTACCTGAGGAGACGGTGACC\" ./origin/*.assembled.fastq |wc")
			
            #os.system("cat ./Igblast_database/human_igh_v ./Igblast_database/human_igh_j ./Igblast_database/human_igh_d > ./Igblast_database/20150429-human-gl-vdj.fasta")
            #os.system("bsub -n 2 -R \"span[hosts=1]\" -q zzh -e err.2 -o out.2 \"2.0.py\"")
            #os.system("bsub -n 1 -q zzh -e err.0 -o out.0 \"misc-histogram-of-seq-length.py\"")
            #os.system("bsub -n 1 -q zzh -e err.0 -o out.0 \"fastqc %s.R1.fastq %s.R2.fastq\""%(file_name, file_name))
            #os.system("bsub -n 8 -q zzh -e err.7 -o out.7 -P %s \"misc_coverage_distribution.py\""%("cover" ))
            
            #os.system(
            #   "bsub -n 8 -q zzh -e err.3 -o out.3 -P %s \"misc_clone_distribution.py\"" %
            #   ("clone"))
            #os.system("bsub -n 8 -q zzh -e err.4 -o out.4 -P %s \"misc_recombanationz_distribution.py\""%("recob" ))
            #os.system("bsub -n 8 -q zzh -e err.5 -o out.5 -P %s \"misc_productive_rate.py\""%("pro" ))
            #os.system("mv ./origin/%s.assembled.fastq  ./origin/%s.fastq"%(file_name, file_name))
            #os.system("bsub -n 2 -q zzh -e err.12 -o out.12 -P %s \"misc_self_balst.py\""%("self_blast" ))
            os.system("bsub -n 2 -q zzh -e err.13 -o out.13 -P %s \"misc_plot_network.py\""%("csvfile" ))
            os.chdir(prj_folder)
    """#deal TCR 
    barcode_sample_dict = {}
    for file_name in os.listdir(prj_folder):
        print file_name
        # and( "K" in file_name or "L" in file_name):
        if os.path.isdir(file_name) :#and "TCR" in file_name:
            os.chdir(file_name)
            num = file_name[-3:]
            print os.getcwd()  # , os.system("which 2.0.py")
            barcode_sample_dict.setdefault(num, []).append(prj_folder + "/" + file_name +"/" + file_name + ".assembled.fastq")
            #os.system("bsub -n 4 -q zzh -e err.11 -o out.11 -P %s \"java -Xmx10G -jar ~/software/mitcr-1.0.3.jar -pset flex -gene TRB ./%s.fastq result_TRB.txt\""%("TCR",file_name))
			
            os.chdir(prj_folder)
    #print barcode_sample_dict
    
    for (key, value) in barcode_sample_dict.items():
        os.system("mkdir sample_TCR_%s"%key)
        all_file_path = " ".join(value)
        #print all_file_path
        os.system("cat %s > ./sample_TCR_%s/sample_TCR_%s.fastq"%(all_file_path, key, key))
    #"""
if __name__ == '__main__':

    pool_size = multiprocessing.cpu_count()
    # create 1st and 2nd subfolders
    prj_folder = os.getcwd()
    prj_tree = ProjectFolders(prj_folder)
    prj_name = fullpath2last_folder(prj_tree.home)
    start = time.time()

    main()
    end = time.time()
    print prj_name, end - start
    print "Finished"
