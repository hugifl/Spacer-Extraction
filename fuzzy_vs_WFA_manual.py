import fuzzysearch
from pywfa import WavefrontAligner
import random
import os
import timeit
import time
from extract_spacers import extract_spacers_wfa
import sys, os, argparse, operator, fuzzysearch
from collections import Counter
import collections
from editSpacer_WFA_old import editSpacer_WFA_old

firstRepeat = 'GAATTGAAAC'
secondRepeat = 'GTCGTACTTT'
fullRepeat = 'GTCGTACTTTACCTAAAAGGAATTGAAAC' #len 29
minStagger = 0
maxStagger = 8
minSpacer = 20
maxSpacer = 60
DR1Mismatch = 2
DR2Mismatch = 2
minReadLength = minStagger + len(firstRepeat) + minSpacer + len(secondRepeat)
# Start timer
dr1_aligner = WavefrontAligner(firstRepeat, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, text_end_free = 29, text_begin_free = 29)
dr2_aligner = WavefrontAligner(secondRepeat, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, text_end_free = 20, text_begin_free = 20)

fasta_directory = "/cluster/scratch/hugifl/fastas/"
files = os.listdir(fasta_directory)

F = open(fasta_directory+files[0], mode='r')
read = F.readlines()[1]

read_short = read[:90]
read = 'AGATGGCCCCTAAAAGGAATTGAAACGCGTTAAGCTCATCAGACAATTTTCAAGCTTATCGGCGTTGACGGGTCGTACTTTACCTAAAAGATTTGTACCAAGGTTCCTAGGACATTACATGATCGGAAGAGCACACGTCTGAACTCCAGTC'
read_1mm = 'AGATGGCCCCTAAAAGGAATTGAAACGCGTTAAGCTCATCAGACAATTTTCAAGCTTATCGGCGTTGACGGGTCGTACTTTACCTAAAAGATTTGTACCAAGGTTCCTAGGACATTACATGATCGGAAGAGCACACGTCTGAACTCCAGTC'+'AAATCCTGTGGGAGCGTGTGGAGCCCTTGAGAGTGC'
read_2mm = 'AGATGGCCCCTAAAAGGAACTGATACGCGTTAAGCTCATCAGACAATTTTCAAGCTTATCGGCGTTGACGGGTCGTACTTTACCTAAAAGATTTGTACCAAGGTTCCTAGGACATTACATGATCGGAAGAGCACACGTCTGAACTCCAGTC'
read_full = 'AGATGGCCCCTAAAAGGAATTGAAACGCGTTAAGCTCATCAGACAATTTTCAAGCTTATCGGCGTTGACGGGTCGTACTTTACCTAAAAGGAATTGAAACAGGTTCCTAGGACATTACATGATCGGAAGAGCACACGTCTGAACTCCAGTC'
test_read = '-------PATTERN1!SSSSSSSSpppppppppaaaaaaaaaaccccceeeeeer!pattern2------------------------------------'
test_read2 = '-------PATTERN1!SSSSSSSSpppppppppaaaaaaaaaaccccceeeeeer!pattern2NNNNNNNNNNNPATTERN1!spaaacccccerrrrererereerr!pattern2------------------------------------'
problem_read = 'AAGGAATTGAAACCAAAATTTTGGGGTTTTCCCCAAAACCCCTTTTGGGGTTTTGTCCGTACTT'
problem_ss1 = 'AAGGAATTGAAACCAAAATTTTGGGGTTTTCCCCAAA'
problem_ss2 = 'CCCCAAAACCCCTTTTGGGGTTTTGTCCGTACTT' 
spacer = "AAAATTTTGGGGTTTTCCCCAAAACCCCTTTTGGGGTTTT"
dr1 = 'GAATTGAAAC'
dr2 = 'GTCGTACTTT'
drf = 'GTCGTACTTTACCTAAAAGGAATTGAAAC'
p1 = 'PATTERN1'
p2 = 'pattern2'
pf = p1+p2
s = 'GGTCATTTATGTCAGACTTGTCGTTTTACAGTTCGTCGTACATTACCTAAAAGATTTGTACCAAGGTTCCTAGATCCACCCAAGATC'
Match = collections.namedtuple('Match', ['start', 'end', 'dist'])


alignment1 = dr1_aligner(problem_ss1, clip_cigar=True, min_aligned_bases_left=2, min_aligned_bases_right=2) #True, min_aligned_bases_left=2, min_aligned_bases_right=2
alignment2 = dr2_aligner(problem_ss2, clip_cigar=True, min_aligned_bases_left=2, min_aligned_bases_right=2)

#print(alignment1.pretty)
#print(alignment.text_start, alignment.text_end)
#print(alignment.aligned_text)
#print(test_read[alignment.text_start + len(p1) -1 : alignment.text_end - len(p2)])
#print(alignment1.text_end)
#print(alignment1.score)
#print(alignment1.text_start, alignment1.text_end)
print(alignment2.pretty)
print(alignment2.score)
#print(alignment2.text_start, alignment2.text_end)
#print(alignment2.text_start, alignment2.text_end)

firstExpect = minStagger
firstRangeToLook = maxStagger - minStagger + len(fullRepeat) - len(firstRepeat)  # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
secondExpect = firstExpect + len(firstRepeat) + minSpacer #30
secondRangeToLook = maxStagger - minStagger + maxSpacer - minSpacer + len(fullRepeat) #58  # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
#spacer=editSpacer_WFA(problem_read,firstRepeat,firstExpect,secondRepeat,secondExpect,dr1_aligner,dr2_aligner,dr2_aligner_dict,minSpacer,maxSpacer,firstRangeToLook,secondRangeToLook)
##
#print(spacer)
#print(len(problem_ss1))
#print(len(problem_ss2))