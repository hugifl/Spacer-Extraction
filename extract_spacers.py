import random
import os
import timeit
from editSpacer_WFA import editSpacer_WFA
from editSpacer_UMI import editSpacer
import sys, os, argparse, operator
from collections import Counter
import collections
import time
from pywfa import WavefrontAligner
import fuzzysearch
from check_double import check_double

# define a custom data type for found matches
Match = collections.namedtuple('Match', ['start', 'end', 'dist'])

def extract_spacers_wfa(F, firstRepeat, secondRepeat, fullRepeat, minReadLength, minSpacer, maxSpacer, dr1_aligner_dict, dr2_aligner_dict, drfull_aligner_dict):
    rawReads = 0
    SinglefullRepeatReads = 0
    spacerReads = 0
    spacerReadsDoubleBoth = 0
    spacerReadsDoubleOne = 0
    full_repeat = 0
    
    iteration_times = []
    start_time_full = time.time()
    for L in F: # loop through reads in file
            if '>' in L: # defline, skip for processing but save read name
                readName = L.strip()
                rawReads+=1
                continue
            
            L=L.strip()
            # ignore reads that are too short to detect adapted spacer
            if len(L) < minReadLength: 
                continue
            
            iteration_start_time = time.time()

            # identify and store reads with more than one acquisition (ie those that contain a full DR sequence)
            skip_read, double, tempSpacerProximal, full_end, L = check_double(L, fullRepeat, drfull_aligner_dict, minReadLength, secondRepeat)
            
            if skip_read:
                continue
    
            if double:
                full_repeat += 1
                spacerProximal = editSpacer_WFA(tempSpacerProximal,firstRepeat,secondRepeat,minSpacer,maxSpacer, dr1_aligner_dict, dr2_aligner_dict)
                #print("proximal: " + spacerProximal) # -------------------------------PRINT ----------------------------------
                #if spacerProximal == '':
                #    print("proximal")
                #    print(tempSpacerProximal)
                #    print("  ")
                tempSpacerDistal = L[(full_end-len(firstRepeat)-8):] # the for is to account for the tendency of regex to add up to 3 nucleotides at the end of a spacer
                
                spacerDistal = ''
                if len(tempSpacerDistal) >= minReadLength:
                    spacerDistal = editSpacer_WFA(tempSpacerDistal,firstRepeat,secondRepeat,minSpacer,maxSpacer, dr1_aligner_dict, dr2_aligner_dict)
                    #print("distal: " + spacerDistal) #-------------------------------PRINT ----------------------------------
                    #if spacerDistal == '':
                    #    print("distal")
                    #    print(tempSpacerDistal)
                    #    print("  ")
                    # if single reads have distal & proximal spacers, label read and export these to file for looking into distal-proximal pairs
                if spacerDistal and spacerProximal:  # !! check how often it is that doubleaq. is found but then 
                    spacerReads+=1
                    spacerReadsDoubleBoth+=1
                # process reads with only distal spacer
                elif spacerDistal:
                    spacerReads+=1
                    spacerReadsDoubleOne+=1
                    

                    
                
                # process reads with only proximal spacer
                elif spacerProximal:
                    spacerReads+=1
                    spacerReadsDoubleOne+=1
                   
                
                else: 
                     
                    continue

            # run standard code if multiple acquisitions (ie full DR sequence) not detected
            else:
                spacer = editSpacer_WFA(L,firstRepeat,secondRepeat,minSpacer,maxSpacer, dr1_aligner_dict, dr2_aligner_dict)
                if spacer == '': 
                    #print(L) # ----------------------------------------print--------------------------------------
                    continue
                spacerReads+=1
                #print(L) # ----------------------------------------print----------------------------------------------
            iteration_end_time = time.time()
            iteration_duration = iteration_end_time - iteration_start_time
            iteration_times.append(iteration_duration) 
                
    end_time_full = time.time()
    total_duration = end_time_full - start_time_full

    return total_duration, iteration_times, rawReads, spacerReads, spacerReadsDoubleBoth, spacerReadsDoubleOne, full_repeat


#def extract_spacers_wfa(F, firstRepeat, secondRepeat, fullRepeat, minReadLength, minSpacer, maxSpacer, dr1_aligner_dict, dr2_aligner_dict, drfull_aligner_dict):
#    rawReads = 0
#    SinglefullRepeatReads = 0
#    spacerReads = 0
#    spacerReadsDoubleBoth = 0
#    spacerReadsDoubleOne = 0
#    full_repeat = 0
#    
#    iteration_times = []
#    start_time_full = time.time()
#    for L in F: # loop through reads in file
#            if '>' in L: # defline, skip for processing but save read name
#                readName = L.strip()
#                rawReads+=1
#                continue
#            
#            L=L.strip()
#            # ignore reads that are too short to detect adapted spacer
#            if len(L) < minReadLength: 
#                continue
#            
#            iteration_start_time = time.time()
#
#            # identify and store reads with more than one acquisition (ie those that contain a full DR sequence)
#            double = False
#            findfull = L.find(fullRepeat)
#            if findfull ==-1:
#                aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
#                drfull_match = aligner_drfull(L, clip_cigar=False)
#                full_start = drfull_match.text_start
#                full_end = drfull_match.text_end
#                double = drfull_match.score > -12 and (drfull_match.text_end - drfull_match.text_start) >= 25
#            elif findfull !=-1:
#                full_start = findfull
#                full_end = findfull + len(fullRepeat)
#                double = True
#         
#            if double:  # if full repeat is present
#                # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
#                tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
#                if len(tempSpacerProximal) < minReadLength:  #          MMAAAAAAAAAAAAAAAABYYYYYYYYY CHAAAAAAAAAANNNNNNNNGGGGGGEEEEEE BBBBBAAAACKKKKKKKK
#                    L = L[len(tempSpacerProximal) - 5:]
#                    double = False
#                    if len(L) >= minReadLength:
#                        aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
#                        drfull_match = aligner_drfull(L, clip_cigar=False)
#                        full_start = drfull_match.text_start
#                        full_end = drfull_match.text_end
#                        double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
#                        tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
#                        if len(tempSpacerProximal) < minReadLength:
#                            double = False
#                    else:
#                        continue
#            
#            if double:
#                full_repeat += 1
#                spacerProximal = editSpacer_WFA(tempSpacerProximal,firstRepeat,secondRepeat,minSpacer,maxSpacer, dr1_aligner_dict, dr2_aligner_dict)
#                #print("proximal: " + spacerProximal) # -------------------------------PRINT ----------------------------------
#                #if spacerProximal == '':
#                #    print("proximal")
#                #    print(tempSpacerProximal)
#                #    print("  ")
#                tempSpacerDistal = L[(full_end-len(firstRepeat)-8):] # the for is to account for the tendency of regex to add up to 3 nucleotides at the end of a spacer
#                
#                spacerDistal = ''
#                if len(tempSpacerDistal) >= minReadLength:
#                    spacerDistal = editSpacer_WFA(tempSpacerDistal,firstRepeat,secondRepeat,minSpacer,maxSpacer, dr1_aligner_dict, dr2_aligner_dict)
#                    #print("distal: " + spacerDistal) #-------------------------------PRINT ----------------------------------
#                    #if spacerDistal == '':
#                    #    print("distal")
#                    #    print(tempSpacerDistal)
#                    #    print("  ")
#                    # if single reads have distal & proximal spacers, label read and export these to file for looking into distal-proximal pairs
#                if spacerDistal and spacerProximal:  # !! check how often it is that doubleaq. is found but then 
#                    spacerReads+=1
#                    spacerReadsDoubleBoth+=1
#                # process reads with only distal spacer
#                elif spacerDistal:
#                    spacerReads+=1
#                    spacerReadsDoubleOne+=1
#                    
#
#                    
#                
#                # process reads with only proximal spacer
#                elif spacerProximal:
#                    spacerReads+=1
#                    spacerReadsDoubleOne+=1
#                    
#                
#                else: 
#                     
#                    continue
#
#            # run standard code if multiple acquisitions (ie full DR sequence) not detected
#            else:
#                spacer = editSpacer_WFA(L,firstRepeat,secondRepeat,minSpacer,maxSpacer, dr1_aligner_dict, dr2_aligner_dict)
#                if spacer == '': 
#                    #print(L) # ----------------------------------------print--------------------------------------
#                    continue
#                spacerReads+=1
#                #print(L) # ----------------------------------------print----------------------------------------------
#            iteration_end_time = time.time()
#            iteration_duration = iteration_end_time - iteration_start_time
#            iteration_times.append(iteration_duration) 
#                
#    end_time_full = time.time()
#    total_duration = end_time_full - start_time_full
#
#    return total_duration, iteration_times, rawReads, spacerReads, spacerReadsDoubleBoth, spacerReadsDoubleOne, full_repeat
#


def extract_spacers_fuzzy(F, firstRepeat, secondRepeat, fullRepeat, minReadLength, minStagger, maxStagger, minSpacer, maxSpacer, DR1Mismatch, DR2Mismatch):
    rawReads = 0
    SinglefullRepeatReads = 0
    spacerReads = 0
    spacerReadsDoubleBoth = 0
    spacerReadsDoubleOne = 0
    full_repeat = 0 
    iteration_times = []
    start_time_full = time.time()

    for L in F: # loop through reads in file
            if '>' in L: # defline, skip for processing but save read name
                readName = L.strip()
                rawReads+=1
                continue
            
            L=L.strip()

            # ignore reads that are too short to detect adapted spacer
            if len(L) < minReadLength: 
                continue

            iteration_start_time = time.time()

            numMismatches = 2
            tempFullRepeat = fuzzysearch.find_near_matches(fullRepeat, L, max_l_dist = numMismatches)
            if tempFullRepeat:  # if full repeat is found

            ## DOUBLE ACQUISITIONS ##
                if len(tempFullRepeat) == 1:  
                    full_repeat += 1
                    SinglefullRepeatReads += 1
                    #I.write(readName+'\n'+L+'\n')

                    # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
                    tempSpacerProximal = L[:tempFullRepeat[0].start+len(secondRepeat)]  
                    firstExpect = minStagger
                    firstRangeToLook = maxStagger - minStagger + len(fullRepeat) - len(firstRepeat) # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                    secondExpect = firstExpect + len(firstRepeat) + minSpacer
                    secondRangeToLook = maxStagger - minStagger + maxSpacer - minSpacer + len(fullRepeat)   # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                    spacerProximal=editSpacer(tempSpacerProximal,firstRepeat,firstExpect,secondRepeat,secondExpect,DR1Mismatch,DR2Mismatch,minSpacer,maxSpacer,firstRangeToLook,secondRangeToLook)
                    #print("proximal: " + spacerProximal) -------------------------------PRINT ----------------------------------

                    tempSpacerDistal = L[tempFullRepeat[0].end-len(firstRepeat):] # the for is to account for the tendency of regex to add up to 3 nucleotides at the end of a spacer
                    firstExpect = 0
                    firstRangeToLook = len(fullRepeat) - len(firstRepeat) # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                    secondExpect = firstExpect + len(firstRepeat) + minSpacer
                    secondRangeToLook = maxSpacer - minSpacer + len(fullRepeat)   # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                    spacerDistal=editSpacer(tempSpacerDistal,firstRepeat,firstExpect,secondRepeat,secondExpect,DR1Mismatch,DR2Mismatch,minSpacer,maxSpacer,firstRangeToLook,secondRangeToLook)
                    #print("distal: " + spacerDistal) --------------------------------------PRINT -------------------------------------
                    # if single reads have distal & proximal spacers, label read and export these to file for looking into distal-proximal pairs

                    if spacerDistal and spacerProximal:
                        spacerReads+=1
                        spacerReadsDoubleBoth+=1

                    
    #
                    # process reads with only distal spacer
                    elif spacerDistal:
                        spacerReads+=1
                        spacerReadsDoubleOne+=1


                    # process reads with only proximal spacer
                    elif spacerProximal:
                        spacerReads+=1
                        spacerReadsDoubleOne+=1

                    else: 
                        continue


            # run standard code if multiple acquisitions (ie full DR sequence) not detected
            else:
                firstExpect = minStagger
                firstRangeToLook = maxStagger - minStagger + len(fullRepeat) - len(firstRepeat)  # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                secondExpect = firstExpect + len(firstRepeat) + minSpacer
                secondRangeToLook = maxStagger - minStagger + maxSpacer - minSpacer + len(fullRepeat)   # this is overkill because using len(fullRepeat), when this could be reduced by knowing the true primer to end of DR distance.
                spacer = editSpacer(L,firstRepeat,firstExpect,secondRepeat,secondExpect,DR1Mismatch,DR2Mismatch,minSpacer,maxSpacer,firstRangeToLook,secondRangeToLook)
                if spacer == '': continue
                spacerReads+=1

            iteration_end_time = time.time()
            iteration_duration = iteration_end_time - iteration_start_time
            iteration_times.append(iteration_duration)

    end_time_full = time.time()
    total_duration = end_time_full - start_time_full

    return total_duration, iteration_times, rawReads, spacerReads, spacerReadsDoubleBoth, spacerReadsDoubleOne, full_repeat
                
    