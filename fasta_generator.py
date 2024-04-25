import random
import sys, os, argparse, operator, fuzzysearch

# Specify the output file name (in the current directory)
output_file = "/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/fake_fastas/fullRepeat_one_spacer.fasta"

spacer = "AAAATTTTGGGGTTTTCCCCAAAACCCCTTTTGGGGTTTT"
firstRepeat = 'GAATTGAAAC'
secondRepeat = 'GTCGTACTTT'
fullRepeat = 'GTCGTACTTTACCTAAAAGGAATTGAAAC'

num_spacer = 0
num_no_spacer = 0
num_dr1 = 0
num_fullR_two_sp = 0
num_fullR_one_sp = 100000


mism_rate = 0.0045
len_empty = 150
len_spacer = 80
len_dr1_only = 150
len_fullR_two_sp = 150
spacer_insert = False
dr1_only_insert = False
spacer_delete = False
dr1_only_delete = False
spacer_var_len = False

def no_spacer_no_dr2_constant_length(firstRepeat, subst_perc, length, insertion, deletion):
    num_nucleotides_before = random.randint(0, 25)
    num_nucleotides_after = length - (num_nucleotides_before + len(firstRepeat))  

    nucleotides_before = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_before))
    nucleotides_after = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_after))

    if insertion:
        insert_pos = random.randint(0, len(firstRepeat))
        firstRepeat = firstRepeat[:insert_pos] + random.choice('ATGC') + firstRepeat[insert_pos:]

    if deletion:
        delete_pos = random.randint(0, len(firstRepeat) - 1)
        firstRepeat = firstRepeat[:delete_pos] +  firstRepeat[delete_pos + 1:]

    if subst_perc > 0:
        firstRepeat = ''.join(
            random.choice('ACGT') if random.random() <= subst_perc else symbol for symbol in firstRepeat)
        
    sequence = f"{nucleotides_before}{firstRepeat}{nucleotides_after}"
    return sequence

def no_spacer_constant_length(length):
    sequence = ''.join(random.choice('ACGT') for _ in range(length))
    return sequence

def single_spacer(spacer, firstRepeat, secondRepeat, subst_perc, length, insertion, deletion, var_len):
    if var_len:
        nu = random.randint(0, 20)
        if random.random() >= 0.5:
            spacer = spacer + spacer[:nu]
        else:
            spacer = spacer[:(len(spacer)-nu)]
            
    
    num_nucleotides_before = random.randint(0, 25)
    num_nucleotides_after = length - (num_nucleotides_before + len(spacer) + len(firstRepeat) + len(secondRepeat))  

    nucleotides_before = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_before))
    nucleotides_after = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_after))

    if insertion:
        insert_pos = random.randint(0, len(firstRepeat))
        if random.random() <= 0.5:
            firstRepeat = firstRepeat[:insert_pos] + random.choice('ATGC') + firstRepeat[insert_pos:]
        else:
            secondRepeat = secondRepeat[:insert_pos] + random.choice('ATGC') + secondRepeat[insert_pos:]

    if deletion:
        delete_pos = random.randint(0, len(firstRepeat) - 1)
        if random.random() <= 0.5:
            firstRepeat = firstRepeat[:delete_pos] +  firstRepeat[delete_pos + 1:]
        else:
            secondRepeat = secondRepeat[:delete_pos] + secondRepeat[delete_pos + 1:]
    
    middle = firstRepeat + spacer + secondRepeat

    if subst_perc > 0:
        middle = ''.join(
            random.choice('ACGT') if random.random() <= subst_perc else symbol for symbol in middle)
        
    sequence = f"{nucleotides_before}{middle}{nucleotides_after}"
    return sequence

def full_repeat_double_spacer_constant_length(spacer, firstRepeat, secondRepeat, fullRepeat, subst_perc, length, insertion, deletion):
    num_nucleotides_before = random.randint(3, 15)
    num_nucleotides_after = length - (num_nucleotides_before + 2 * len(spacer) + len(firstRepeat) + len(secondRepeat) + len(fullRepeat))  

    nucleotides_before = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_before))
    nucleotides_after = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_after))

    if insertion:
        insert_pos = random.randint(0, len(firstRepeat))
        insert_pos_2 = random.randint(0, len(fullRepeat))
        if random.random() <= 0.5:
            firstRepeat = firstRepeat[:insert_pos] + random.choice('ATGC') + firstRepeat[insert_pos:]
        else:
            secondRepeat = secondRepeat[:insert_pos] + random.choice('ATGC') + secondRepeat[insert_pos:]
        fullRepeat = fullRepeat[:insert_pos_2] + random.choice('ATGC') + fullRepeat[insert_pos_2:]

    if deletion:
        delete_pos = random.randint(0, len(firstRepeat) - 1)
        delete_pos_2 = random.randint(0, len(fullRepeat) - 1)
        if random.random() <= 0.5:
            firstRepeat = firstRepeat[:delete_pos] +  firstRepeat[delete_pos + 1:]
        else:
            secondRepeat = secondRepeat[:delete_pos] + secondRepeat[delete_pos + 1:]
        fullRepeat = fullRepeat[:delete_pos_2] + fullRepeat[delete_pos_2 + 1:]
    

    if subst_perc > 0:
        fullRepeat = ''.join(
            random.choice('ACGT') if random.random() <= subst_perc else symbol for symbol in fullRepeat)
        firstRepeat = ''.join(
            random.choice('ACGT') if random.random() <= subst_perc else symbol for symbol in firstRepeat)
        secondRepeat = ''.join(
            random.choice('ACGT') if random.random() <= subst_perc else symbol for symbol in secondRepeat)
        
    sequence = f"{nucleotides_before}{firstRepeat}{spacer}{fullRepeat}{spacer}{secondRepeat}{nucleotides_after}"
    return sequence

def full_repeat_single_spacer_constant_length(spacer, firstRepeat, secondRepeat, fullRepeat, subst_perc, length, insertion, deletion):
    num_nucleotides_before = random.randint(3, 15)
    num_nucleotides_after = length - (num_nucleotides_before + len(spacer) + len(firstRepeat) + len(fullRepeat))  

    nucleotides_before = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_before))
    nucleotides_after = ''.join(random.choice('ACGT') for _ in range(num_nucleotides_after))

    if insertion:
        insert_pos = random.randint(0, len(firstRepeat))
        insert_pos_2 = random.randint(0, len(fullRepeat))
        if random.random() <= 0.5:
            firstRepeat = firstRepeat[:insert_pos] + random.choice('ATGC') + firstRepeat[insert_pos:]
        else:
            secondRepeat = secondRepeat[:insert_pos] + random.choice('ATGC') + secondRepeat[insert_pos:]
        fullRepeat = fullRepeat[:insert_pos_2] + random.choice('ATGC') + fullRepeat[insert_pos_2:]

    if deletion:
        delete_pos = random.randint(0, len(firstRepeat) - 1)
        delete_pos_2 = random.randint(0, len(fullRepeat) - 1)
        if random.random() <= 0.5:
            firstRepeat = firstRepeat[:delete_pos] +  firstRepeat[delete_pos + 1:]
        else:
            secondRepeat = secondRepeat[:delete_pos] + secondRepeat[delete_pos + 1:]
        fullRepeat = fullRepeat[:delete_pos_2] + fullRepeat[delete_pos_2 + 1:]
    

    if subst_perc > 0:
        fullRepeat = ''.join(
            random.choice('ACGT') if random.random() <= subst_perc else symbol for symbol in fullRepeat)
        firstRepeat = ''.join(
            random.choice('ACGT') if random.random() <= subst_perc else symbol for symbol in firstRepeat)
        secondRepeat = ''.join(
            random.choice('ACGT') if random.random() <= subst_perc else symbol for symbol in secondRepeat)
        
    sequence = f"{nucleotides_before}{firstRepeat}{spacer}{fullRepeat}{nucleotides_after}"
    return sequence

spacer_seq = [single_spacer(spacer, firstRepeat, secondRepeat, mism_rate, len_spacer, spacer_insert, spacer_delete, spacer_var_len) for _ in range(num_spacer)]
empty_seq = [no_spacer_constant_length(len_empty) for _ in range(num_no_spacer)]
only_dr1_seq = [no_spacer_no_dr2_constant_length(firstRepeat, mism_rate, len_dr1_only, dr1_only_insert, dr1_only_delete) for _ in range(num_dr1)]
fullRep_two_spacer = [full_repeat_double_spacer_constant_length(spacer, firstRepeat, secondRepeat,fullRepeat, mism_rate, len_fullR_two_sp, spacer_insert, spacer_delete) for _ in range(num_fullR_two_sp)]
fullRep_one_spacer = [full_repeat_single_spacer_constant_length(spacer, firstRepeat, secondRepeat,fullRepeat, mism_rate, len_fullR_two_sp, spacer_insert, spacer_delete) for _ in range(num_fullR_one_sp)]


sequences = fullRep_two_spacer + spacer_seq + empty_seq + only_dr1_seq + fullRep_one_spacer


# Write the sequences to the FASTA file in the current directory

with open(output_file, "w") as fasta_file:
    for sequence in sequences:
        fasta_file.write(f"{sequence}\n")

print("done")