import fuzzysearch
from pywfa import WavefrontAligner
import random
import os
import timeit
import time
#from extract_spacers import extract_spacers_wfa
import sys, os, argparse, operator, fuzzysearch
from collections import Counter
import collections
from WFA2_aligner_construction import construct_aligners


firstRepeat = 'GAATTGAAAC'
test_seq = 'TESTSEQUENCE'
secondRepeat = 'GTCGTACTTT'
both_repeat = firstRepeat + secondRepeat
spacer = 'AAAATTTTGGGGAAAACCCCTTTTAAAACCCCTTTTGGGG'
fullRepeat = 'GTCGTACTTTACCTAAAAGGAATTGAAAC' #len 29
minStagger = 0
maxStagger = 8
minSpacer = 20
maxSpacer = 60
DR1Mismatch = 2
DR2Mismatch = 2
minReadLength = minStagger + len(firstRepeat) + minSpacer + len(secondRepeat)

#fuzzysearch
def fuzzy_align(pattern, text, max_dist):
    res = fuzzysearch.find_near_matches(pattern, text, max_l_dist = max_dist)
    res = sorted(res, key=lambda x: x.dist)
    return res[0].start, res[0].end
##pywfa

# construct the wfa2 alignment objects 
dr1_aligner_dict = construct_aligners(firstRepeat, 10, 150)  
dr2_aligner_dict = construct_aligners(secondRepeat, 10, 150) 
drfull_aligner_dict = construct_aligners(fullRepeat, 10, 150)
dr1_aligner_error = WavefrontAligner(firstRepeat, gap_opening = 4, gap_opening2 = 24, gap_extension = 4, text_end_free = 0, text_begin_free = 0)
#dr_full_aligner = WavefrontAligner(drf)
#pattern_aligner = WavefrontAligner(pf, gap_opening = 0, gap_opening2 = 600, mismatch = 400, gap_extension = 1)
#complete_aligner = WavefrontAligner(dr1+dr2, gap_opening = 0, gap_opening2 = 150000, gap_extension = 1, mismatch = 100000)

run_length = False
run_length_error = False
run_pos = False
run_dr1_spacer_dr2 = False
run_length_find = False
run_spacer_alignment = False
run_spacer_alignment_match = False
run_spacer_alignment_match_2 = True
run_spacer_alignment_match_subspacer_length = False

def dr1_sequence(sequence, MM_no, length, position):
    start_position = int(((position * length)-(0.5 * len(sequence))))     
    num_nucleotides_before = start_position
    num_nucleotides_after = length - (num_nucleotides_before + len(sequence))  

    nucleotides_before = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_before))
    nucleotides_after = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_after))

    MM_positions = set()
    for _ in range(MM_no):
        while True:
            position = random.randint(0, len(sequence) - 1)
            if position not in MM_positions:
                MM_positions.add(position)
                break

        letter = random.choice('ATGC')
        sequence = sequence[:position] + letter + sequence[position + 1:]

    read = f"{nucleotides_before}{sequence}{nucleotides_after}"
    return read

def dr1_spacer_dr2sequence(dr1, spacer, dr2, MM_no, length, position):
    MM_positions = set()
    for _ in range(MM_no):
        while True:
            MM_position = random.randint(0, len(dr1) - 1)
            if MM_position not in MM_positions:
                MM_positions.add(MM_position)
                break

        letter = random.choice('ATGC')
        dr1 = dr1[:MM_position] + letter + dr1[MM_position + 1:]
    
    MM_positions = set()
    for _ in range(MM_no):
        while True:
            MM_position = random.randint(0, len(dr2) - 1)
            if MM_position not in MM_positions:
                MM_positions.add(MM_position)
                break

        letter = random.choice('ATGC')
        dr2 = dr2[:MM_position] + letter + dr2[MM_position + 1:]

    middle = dr1 + spacer + dr2
    
    start_position = int(((position * length)-(0.5 * len(middle))))     
    num_nucleotides_before = 5
    num_nucleotides_after = length - (num_nucleotides_before + len(middle))  
    nucleotides_before = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_before))
    nucleotides_after = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_after))

    read = f"{nucleotides_before}{middle}{nucleotides_after}"
    return read

def matching_spacers(MM_no, length):
    
    spacer_text = ''.join(random.choice('ACGT') for _ in range(length))
    spacer_pattern = spacer_text
    MM_positions = set()
    for _ in range(MM_no):
        while True:
            position = random.randint(0, len(spacer_pattern) - 1)
            if position not in MM_positions:
                MM_positions.add(position)
                break

        letter = random.choice('ATGC')
        spacer_pattern = spacer_pattern[:position] + letter + spacer_pattern[position + 1:]
    return spacer_text, spacer_pattern

def matching_spacers_partial(MM_no, length):
    spacer_text = ''.join(random.choice('ACGT') for _ in range(120))
    spacer_pattern = spacer_text
    MM_positions = set()
    for _ in range(MM_no):
        while True:
            position = random.randint(0, len(spacer_pattern) - 1)
            if position not in MM_positions:
                MM_positions.add(position)
                break

        letter = random.choice('ATGC')
        spacer_pattern = spacer_pattern[:position] + letter + spacer_pattern[position + 1:]
    start = (60-int(0.5*length))
    end = (60+int(0.5*length))
    spacer_pattern = spacer_pattern[start:end]
    return spacer_text, spacer_pattern

def random_spacers(length):
    nucleotides = ''.join(random.choice('ACGT') for _ in range(length))
    return nucleotides

if run_spacer_alignment_match_subspacer_length:
    Method = []
    Length = []
    Time = []
    l_dist = []
    No_MM = []

    for MM_no in range(6):
        print(MM_no)
        for length in range(20,121): #121
            print(length)
            spacers = [matching_spacers_partial(MM_no, length) for _ in range(100)] 
            #fuzzy align score = 2
            start_time = time.time()
            for spacer_pair in spacers:
                spacer_text = spacer_pair[0]
                spacer_pattern = spacer_pair[1]
                match = fuzzysearch.find_near_matches(spacer_pattern, spacer_text, max_l_dist = 2) 
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            Method.append('fuzzy')
            Length.append(length)
            l_dist.append("max_l_dist = 2")
            No_MM.append(MM_no)
            #fuzzy align score = 4
            start_time = time.time()
            for spacer_pair in spacers:
                spacer_text = spacer_pair[0]
                spacer_pattern = spacer_pair[1]
                match = fuzzysearch.find_near_matches(spacer_pattern, spacer_text, max_l_dist = 4)
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            Method.append('fuzzy')
            Length.append(length)
            l_dist.append('max_l_dist = 4')
            No_MM.append(MM_no)

            #fuzzy align score = 6
            start_time_a = time.time()
            for spacer_pair in spacers:
                spacer_text = spacer_pair[0]
                spacer_pattern = spacer_pair[1]
                match = fuzzysearch.find_near_matches(spacer_pattern, spacer_text, max_l_dist = 6)
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            Method.append('fuzzy')
            Length.append(length)
            l_dist.append("max_l_dist = 6")
            No_MM.append(MM_no)
            

    csv_file = "/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/time_experiments_outputs/spacer_alignment_match_partial.csv"
    data = list(zip(Method, l_dist, Length, No_MM, Time))
    with open(csv_file, 'w', newline='') as file:
        import csv
        writer = csv.writer(file)
        writer.writerow(["Method", "score_cutoff", "Length", "No_MM", "Time"]) 
        writer.writerows(data)
    print(f"CSV file '{csv_file}' created.")





if run_spacer_alignment_match:
    Method = []
    Length = []
    Time = []
    l_dist = []
    No_MM = []

    for MM_no in range(6):
        for length in range(20,121): #121
            spacers = [matching_spacers(MM_no, length) for _ in range(100)] 
            #fuzzy align score = 2
            start_time = time.time()
            for spacer_pair in spacers:
                spacer_text = spacer_pair[0]
                spacer_pattern = spacer_pair[1]
                match = fuzzysearch.find_near_matches(spacer_pattern, spacer_text, max_l_dist = 2) 
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            Method.append('fuzzy')
            Length.append(length)
            l_dist.append("max_l_dist = 2")
            No_MM.append(MM_no)
            #fuzzy align score = 4
            start_time = time.time()
            for spacer_pair in spacers:
                spacer_text = spacer_pair[0]
                spacer_pattern = spacer_pair[1]
                match = fuzzysearch.find_near_matches(spacer_pattern, spacer_text, max_l_dist = 4)
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            Method.append('fuzzy')
            Length.append(length)
            l_dist.append('max_l_dist = 4')
            No_MM.append(MM_no)

            #fuzzy align score = 6
            start_time_a = time.time()
            for spacer_pair in spacers:
                spacer_text = spacer_pair[0]
                spacer_pattern = spacer_pair[1]
                match = fuzzysearch.find_near_matches(spacer_pattern, spacer_text, max_l_dist = 6)
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            Method.append('fuzzy')
            Length.append(length)
            l_dist.append("max_l_dist = 6")
            No_MM.append(MM_no)
            #WFA align
            start_time = time.time()
            for spacer_pair in spacers:
                spacer_text = spacer_pair[0]
                spacer_pattern = spacer_pair[1]
                aligner = WavefrontAligner(spacer_pattern, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                match = aligner(spacer_text, clip_cigar=False)  
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            Method.append('wfa')
            Length.append(length)
            l_dist.append("-")
            No_MM.append(MM_no)

    csv_file = "/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/time_experiments_outputs/spacer_alignment_match.csv"
    data = list(zip(Method, l_dist, Length, No_MM, Time))
    with open(csv_file, 'w', newline='') as file:
        import csv
        writer = csv.writer(file)
        writer.writerow(["Method", "score_cutoff", "Length", "No_MM", "Time"]) 
        writer.writerows(data)
    print(f"CSV file '{csv_file}' created.")

if run_spacer_alignment:
    Method = []
    Length = []
    Time = []
    l_dist = []
    spacer = 'AAAACCCCTTTTGGGGAAAAAAAACCCCTTTTGGGGAAAA' 
    
    for length in range(20,121):
        spacers = [random_spacers(length) for _ in range(1000)] 
        #fuzzy align score = 2
        start_time = time.time()
        for L in spacers:
            match = fuzzysearch.find_near_matches(spacer, L, max_l_dist = 2) 
        end_time = time.time()
        measure_time = end_time - start_time
        Time.append(measure_time)
        Method.append('fuzzy')
        Length.append(length)
        l_dist.append("max_l_dist = 2")
        #fuzzy align score = 4
        start_time = time.time()
        for L in spacers:
            match = fuzzysearch.find_near_matches(spacer, L, max_l_dist = 4) 
        end_time = time.time()
        measure_time = end_time - start_time
        Time.append(measure_time)
        Method.append('fuzzy')
        Length.append(length)
        l_dist.append('max_l_dist = 4')
        
        #fuzzy align score = 6
        start_time = time.time()
        for L in spacers:
            match = fuzzysearch.find_near_matches(spacer, L, max_l_dist = 6) 
        end_time = time.time()
        measure_time = end_time - start_time
        Time.append(measure_time)
        Method.append('fuzzy')
        Length.append(length)
        l_dist.append("max_l_dist = 6")
        #WFA align
        start_time = time.time()
        for L in spacers:
            aligner = WavefrontAligner(spacer, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
            match = aligner(L, clip_cigar=False)  
        end_time = time.time()
        measure_time = end_time - start_time
        Time.append(measure_time)
        Method.append('wfa')
        Length.append(length)
        l_dist.append("-")

    csv_file = "/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/time_experiments_outputs/spacer_alignment.csv"
    data = list(zip(Method, l_dist, Length, Time))
    with open(csv_file, 'w', newline='') as file:
        import csv
        writer = csv.writer(file)
        writer.writerow(["Method", "score_cutoff", "Length","Time"]) 
        writer.writerows(data)
    print(f"CSV file '{csv_file}' created.")

if run_dr1_spacer_dr2:
    No_MM = []
    Method = []
    Length = []
    Time = []

    for MM in range(1,3):
        for length in range(75,151):
            sequences = [dr1_spacer_dr2sequence(firstRepeat,spacer, secondRepeat,MM,length,0.5) for _ in range(1000)] 
            #fuzzy align
            start_time = time.time()
            for L in sequences:
                Dr1_matches = fuzzysearch.find_near_matches(firstRepeat, L[:30], max_l_dist = 2) #30 nu long is the length of the string where dr2 is searched in the current code for a 150 nu read.
                Dr2_matches = fuzzysearch.find_near_matches(secondRepeat, L[5:73], max_l_dist = 2) #68 nu long is the length of the string where dr2 is searched in the current code for a 150 nu read.
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            No_MM.append(MM)
            Method.append('fuzzy')
            Length.append(length)

            #WFA align
            start_time = time.time()
            for L in sequences:
                dr1_aligner = dr1_aligner_dict[round(len(L), -1)-10]
                dr1_match = dr1_aligner(L, clip_cigar=False)
                dr2_aligner = dr2_aligner_dict[round(len(L), -1)-10]
                dr2_match = dr2_aligner(L, clip_cigar=False)
                    
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            No_MM.append(MM)
            Method.append('wfa')
            Length.append(length)

    csv_file = "/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/time_experiments_outputs/length_dr1_spacer_dr2.csv"
    data = list(zip(Method, No_MM, Length, Time))
    with open(csv_file, 'w', newline='') as file:
        import csv
        writer = csv.writer(file)
        writer.writerow(["Method", "No_MM", "Length","Time"]) 
        writer.writerows(data)
    print(f"CSV file '{csv_file}' created.")

if run_length_find:
   
    Method = []
    Length = []
    Time = []

    
    for length in range(20,151):
        sequences = [dr1_sequence(firstRepeat,0,length,0.5) for _ in range(1000)] 
        #string.find
        start_time = time.time()
        for L in sequences:
            findDR1 = L.find(firstRepeat)
        end_time = time.time()
        measure_time = end_time - start_time
        Time.append(measure_time)
        Method.append('string_find')
        Length.append(length)
        #WFA align
        start_time = time.time()
        for L in sequences:
            aligner = dr1_aligner_dict[round(len(L), -1)-10]
            Match = aligner(L, clip_cigar=False)
        end_time = time.time()
        measure_time = end_time - start_time
        Time.append(measure_time)
        Method.append('wfa')
        Length.append(length)

    csv_file = "/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/time_experiments_outputs/length_find.csv"
    data = list(zip(Method, Length, Time))
    with open(csv_file, 'w', newline='') as file:
        import csv
        writer = csv.writer(file)
        writer.writerow(["Method", "Length","Time"]) 
        writer.writerows(data)
    print(f"CSV file '{csv_file}' created.")


if run_length:
    No_MM = []
    Method = []
    Length = []
    Time = []

    for MM in range(1,3):
        for length in range(20,151):
            sequences = [dr1_sequence(firstRepeat,MM,length,0.5) for _ in range(1000)] 

            #fuzzy align
            start_time = time.time()
            for L in sequences:
                Matches = fuzzysearch.find_near_matches(firstRepeat, L, max_l_dist = 2)
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            No_MM.append(MM)
            Method.append('fuzzy')
            Length.append(length)

            #WFA align
            start_time = time.time()
            for L in sequences:
                aligner = dr1_aligner_dict[round(len(L), -1)-10]
                Match = aligner(L, clip_cigar=False)
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            No_MM.append(MM)
            Method.append('wfa')
            Length.append(length)

    #csv_file = "/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/time_experiments_outputs/length.csv"
    #data = list(zip(Method, No_MM, Length, Time))
    #with open(csv_file, 'w', newline='') as file:
    #    import csv
    #    writer = csv.writer(file)
    #    writer.writerow(["Method", "No_MM", "Length","Time"]) 
    #    writer.writerows(data)
    #print(f"CSV file '{csv_file}' created.")

if run_spacer_alignment_match_2:
    Method = []
    Length = []
    Time = []
    l_dist = []
    No_MM = []

    for MM_no in [0, 2, 4]:
        print("MM: ",MM_no)
        for length in range(50,161): #121
            sequences = [dr1_sequence(firstRepeat,MM_no,length,0.5) for _ in range(1000)]  
            #fuzzy align score = 2
            start_time = time.time()
            for L in sequences:
                Matches = fuzzysearch.find_near_matches(firstRepeat, L, max_l_dist = 2) 
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            Method.append('fuzzy')
            Length.append(length)
            l_dist.append(f"max_l_dist = {2}")
            No_MM.append(MM_no)
            #fuzzy align score = 4
            start_time = time.time()
            for L in sequences:
                Matches = fuzzysearch.find_near_matches(firstRepeat, L, max_l_dist = 4) 
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            Method.append('fuzzy')
            Length.append(length)
            l_dist.append('max_l_dist = 4')
            No_MM.append(MM_no)

            #fuzzy align score = 6
            start_time_a = time.time()
            for L in sequences:
                Matches = fuzzysearch.find_near_matches(firstRepeat, L, max_l_dist = 6) 
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            Method.append('fuzzy')
            Length.append(length)
            l_dist.append("max_l_dist = 6")
            No_MM.append(MM_no)
            #WFA align correct
            start_time = time.time()
            for L in sequences:
                aligner = dr1_aligner_dict[round(len(L), -1)-10]
                Match = aligner(L, clip_cigar=False) 
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            Method.append('wfa')
            Length.append(length)
            l_dist.append("-")
            No_MM.append(MM_no)

            #WFA align wrong
            start_time = time.time()
            for L in sequences:
                Match = dr1_aligner_error(L, clip_cigar=False)
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            Method.append('wfa')
            Length.append(length)
            l_dist.append("--")
            No_MM.append(MM_no)


    csv_file = "/cluster/home/hugifl/dev-hugi/time_experiments_outputs/spacer_alignment_length_wfa_dict_vs_error.csv"
    data = list(zip(Method, l_dist, Length, No_MM, Time))
    with open(csv_file, 'w', newline='') as file:
        import csv
        writer = csv.writer(file)
        writer.writerow(["Method", "score_cutoff", "Length", "No_MM", "Time"]) 
        writer.writerows(data)
    print(f"CSV file '{csv_file}' created.")


if run_length_error:
    No_MM = []
    Method = []
    Length = []
    Time = []

    for MM in range(1,3):
        for length in range(20,151):
            sequences = [dr1_sequence(firstRepeat,MM,length,0.5) for _ in range(1000)] 

            #fuzzy align
            start_time = time.time()
            for L in sequences:
                Matches = fuzzysearch.find_near_matches(firstRepeat, L, max_l_dist = 2)
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            No_MM.append(MM)
            Method.append('fuzzy')
            Length.append(length)

            #WFA align
            start_time = time.time()
            for L in sequences:
                Match = dr1_aligner_error(L, clip_cigar=False)
            end_time = time.time()
            measure_time = end_time - start_time
            Time.append(measure_time)
            No_MM.append(MM)
            Method.append('wfa')
            Length.append(length)

    csv_file = "/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/time_experiments_outputs/length_error.csv"
    data = list(zip(Method, No_MM, Length, Time))
    with open(csv_file, 'w', newline='') as file:
        import csv
        writer = csv.writer(file)
        writer.writerow(["Method", "No_MM", "Length","Time"]) 
        writer.writerows(data)
    print(f"CSV file '{csv_file}' created.")



if run_pos:
    No_MM = []
    Method = []
    Length = []
    Time = []
    Position = []

    for MM in range(1,2): 
        for pos in range(0,11):
            for length in range(30,151,20): 
                sequences = [dr1_sequence(firstRepeat,MM,length,pos/10) for _ in range(1000)] 

                #fuzzy align
                start_time = time.time()
                for L in sequences:
                    Matches = fuzzysearch.find_near_matches(firstRepeat, L, max_l_dist = 2)
                end_time = time.time()
                measure_time = end_time - start_time
                Time.append(measure_time)
                No_MM.append(MM)
                Method.append('fuzzy')
                Length.append(length)
                Position.append(pos)
                #WFA align
                start_time = time.time()
                for L in sequences:
                    aligner = dr1_aligner_dict[round(len(L), -1)-10]
                    Match = aligner(L, clip_cigar=False)
                end_time = time.time()
                measure_time = end_time - start_time
                Time.append(measure_time)
                No_MM.append(MM)
                Method.append('wfa')
                Length.append(length)
                Position.append(pos)

    csv_file = "/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/time_experiments_outputs/position.csv"
    data = list(zip(Method, No_MM, Length, Position, Time))
    with open(csv_file, 'w', newline='') as file:
        import csv
        writer = csv.writer(file)
        writer.writerow(["Method", "No_MM", "Length","Position","Time"]) 
        writer.writerows(data)
    print(f"CSV file '{csv_file}' created.")