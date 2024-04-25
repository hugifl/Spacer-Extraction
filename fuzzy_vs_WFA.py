import fuzzysearch
from pywfa import WavefrontAligner
import random
import numpy as np
import os
import timeit
from editSpacer_UMI import editSpacer
import sys, argparse, operator, fuzzysearch
from collections import Counter
import collections
import time
from extract_spacers import extract_spacers_fuzzy, extract_spacers_wfa
from WFA_aligner_construction import construct_aligners
import gzip

fasta_directory = "/cluster/scratch/hugifl/paraquat_run_4/intermediates/fastas/" 

files = os.listdir(fasta_directory)

#fake_fasta_directory = "/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/fake_fastas/"
#fake_files = os.listdir(fake_fasta_directory) 
fasta_file_path = fasta_directory+files[0]


F = open(fasta_file_path, mode='r')
#F = open(fake_fasta_directory+fake_files[10], mode='r')
#print(files[0])

            

# define a custom data type for found matches
Match = collections.namedtuple('Match', ['start', 'end', 'dist'])

UMI_status ='off'


UMI_length = 10
#LBC = 
firstRepeat = 'GAATTGAAAC'
secondRepeat = 'GTCGTACTTT'
fullRepeat = 'GTCGTACTTTACCTAAAAGGAATTGAAAC' #len 29
LBC = 'TTTGTACCAAGGTTCCTAG'
minStagger = 0
maxStagger = 8
minSpacer = 20
maxSpacer = 60
DR1Mismatch = 2
DR2Mismatch = 0
minReadLength = minStagger + len(firstRepeat) + minSpacer + len(secondRepeat) #40


#def construct_aligners(pattern, min_text_length, max_text_length):
#    aligner_dict = {}
#    for length in range(min_text_length, max_text_length + 1, 10):
#        aligner = WavefrontAligner(pattern, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, text_end_free = length, text_begin_free = length)
#        aligner_dict[length] = aligner
#
#    return aligner_dict
#
#dr1_aligner_dict = construct_aligners(firstRepeat, 20, 150)  #dr1_aligners()
#dr2_aligner_dict = construct_aligners(secondRepeat, 20, 150) #dr2_aligners()
#drfull_aligner_dict = construct_aligners(fullRepeat, 20, 150)#drfull_aligners()

# construct the wfa2 alignment objects 
dr1_aligner_dict = construct_aligners(firstRepeat, 20, 150)  
dr2_aligner_dict = construct_aligners(secondRepeat, 20, 150) 
drfull_aligner_dict = construct_aligners(fullRepeat, 20, 150)


fuzzy =   False
wfa =  True

if fuzzy:
    t_tot, t_its, rawReads, spacerReads, spacerReadsDoubleBoth, spacerReadsDoubleOne, full_repeat = extract_spacers_fuzzy(F, firstRepeat, secondRepeat, fullRepeat, minReadLength, minStagger, maxStagger, minSpacer, maxSpacer, DR1Mismatch, DR2Mismatch)
    mean = np.mean(t_its)
    sd = np.std(t_its)
    print("---------------------- FUZZY ALIGNER ----------------------")
    print(f"Total time: {t_tot}\nAverage time per read: {mean}\nStd of time per read: {sd}\nRaw reads: {rawReads}\nSpacerReads: {spacerReads}\nSpacer reads double both: {spacerReadsDoubleBoth}\nSpacer reads double one: {spacerReadsDoubleOne}\nFull_repeats : {full_repeat}")


if wfa:
    t_tot_3, t_its_3, rawReads_3, spacerReads_3, spacerReadsDoubleBoth_3, spacerReadsDoubleOne_3, full_repeat_3 = extract_spacers_wfa(F, firstRepeat, secondRepeat, fullRepeat, minReadLength, minSpacer, maxSpacer, dr1_aligner_dict, dr2_aligner_dict, drfull_aligner_dict)
    mean_3 = np.mean(t_its_3)
    sd_3 = np.std(t_its_3)
    print("---------------------- WFA ALIGNER ----------------------")
    print(f"Total time: {t_tot_3}\nAverage time per read: {mean_3}\nStd of time per read: {sd_3}\nRaw reads: {rawReads_3}\nSpacerReads: {spacerReads_3}\nSpacer reads double both: {spacerReadsDoubleBoth_3}\nSpacer reads double one: {spacerReadsDoubleOne_3}\nFull_repeats : {full_repeat_3}")
