import matplotlib
from pywfa import WavefrontAligner
import math, random, os
import numpy as np
import collections, fuzzysearch
from WFA2_aligner_construction import construct_aligners
from editSpacer_WFA import editSpacer_WFA
import csv, sys

#################### Which experiment to run #####################
experiment_1 = False
experiment_2 = False
experiment_3 = False
experiment_4 = False
experiment_5 = False
experiment_6 = False
experiment_6_5 = True
experiment_7 = False
experiment_8 = False
experiment_8_5 = False
########### Which fasta file to run the experiment on ############
fasta_file = 5
#fasta_files = [0] #[0,5,10,15,20,25,30]




fasta_directory = "/cluster/scratch/hugifl/fastas/"
files = os.listdir(fasta_directory)
fake_fasta_directory = "/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/fake_fastas/"
fake_files = os.listdir(fake_fasta_directory) 


firstRepeat = 'GAATTGAAAC'
secondRepeat = 'GTCGTACTTT'
secondRepeat_full = 'GTCGTACTTTACCTAAAAG'
fullRepeat = 'GTCGTACTTTACCTAAAAGGAATTGAAAC'
LBC = 'ATTTGTACCAAGGTTCCTAG'

# construct the wfa2 alignment objects 
dr1_aligner_dict = construct_aligners(firstRepeat, 10, 150)  
dr2_aligner_dict = construct_aligners(secondRepeat, 10, 150) 
drfull_aligner_dict = construct_aligners(fullRepeat, 10, 150)
LBC_aligner_dict = construct_aligners(LBC, 20, 150)


Match = collections.namedtuple('Match', ['start', 'end', 'dist'])

F = open(fasta_directory+files[fasta_file], mode='r')




################## EXPERIMENT 1: Here I test how many of the identified double acquisitions (full read present) have a LBC using WFA2 alignments #########################
#
# Findings: None of the real double acquisitions have a LBC since the reads are too short. All cases where a LBC is found are single acquisitions
#           where the Dr2_LBC is close enough to a full repeat (Dr2Dr1) so that it gets identified as a double acquisition.
#           If a stricter score cut off (>-14) is used to find full repeats, 8885 full repeats are identified, 8830 of which have no LBC. The 55 which have a LBC
#           are FP full repeats. If a less stirct score cut off (>-16) is used, ca. 14000 full repeats are found, 9000 of which have no LBC. The ca. 5000 that have a LBC
#           are FP full repeats.
if experiment_1:
    txt_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_1_{files[fasta_file]}"
    not_found = 0
    found = 0
    tot = 0
    full_repeat = 0
    for L in F:
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            continue
        
        L=L.strip()
        # ignore reads that are too short to detect adapted spacer
        if len(L) < 60: 
            continue
        tot += 1

        double = False
        findfull = L.find(fullRepeat)
        if findfull ==-1:
            aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
            drfull_match = aligner_drfull(L, clip_cigar=False)
            full_start = drfull_match.text_start
            full_end = drfull_match.text_end
            double = drfull_match.score > -16 and (drfull_match.text_end - drfull_match.text_start) >= 25
        elif findfull !=-1:
            full_start = findfull
            full_end = findfull + len(fullRepeat)
            double = True

        if double:  # if full repeat is present
            # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
            tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
            if len(tempSpacerProximal) < 60:  
                L = L[len(tempSpacerProximal) - 5:]
                double = False
                if len(L) >= 60:
                    aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                    drfull_match = aligner_drfull(L, clip_cigar=False)
                    full_start = drfull_match.text_start
                    full_end = drfull_match.text_end
                    double = drfull_match.score > -16 and (drfull_match.text_end - drfull_match.text_start) >= 25
                    tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
                    if len(tempSpacerProximal) < 60:
                        double = False
                else:
                    continue
                
        if double:
            full_repeat += 1
            LBC_end = ''
            findLBC = L.find(LBC)
            if findLBC==-1:
                aligner_LBC = LBC_aligner_dict[(len(L) // 10) * 10]
                LBC_match = aligner_LBC(L, clip_cigar=False)
                if (LBC_match.text_end - LBC_match.text_start)>= 6 and LBC_match.score > -14:             
                    LBC_end = LBC_match.text_end




            else:
                LBC_match = [Match(findLBC, findLBC+len(LBC),0)]
                LBC_end = LBC_match[0].end    
            if not LBC_end:
                not_found += 1
            if LBC_end:
                print(drfull_match.score)
                print(drfull_match.pretty)





        else:
            continue

    with open(txt_file, 'w') as file:
        sys.stdout = file
        print("tot reads: "+str(tot))
        print("full repeats without LBC: "+str(not_found))
        print("full repeats totally: "+str(full_repeat))
    sys.stdout = sys.__stdout__
    print("tot reads: "+str(tot))
    print("full repeats without LBC: "+str(not_found))
    print("full repeats totally: "+str(full_repeat))


################## EXPERIMENT 2: Here I test how many of the reads that have a LBC are double acquisitions (fullrepeat present) using WFA2 alignments #########################
#
# Findings: (Amost) none of the full repeats identified after restricting for the presence of LBC have a real full repeat. All cases where a LBC is found are single acquisitions
#           where the Dr2_LBC is close enough to a full repeat (Dr2Dr1) so that it gets identified as a double acquisition.
#
if experiment_2:
    txt_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_2_{files[fasta_file]}"
    not_found = 0
    tot = 0
    full_repeat = 0
    full_early = 0
    for L in F:
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            continue
        
        L=L.strip()
        # ignore reads that are too short to detect adapted spacer
        if len(L) < 60: 
            continue
        tot += 1
        LBC_end = ''
        findLBC = L.find(LBC)
        if findLBC==-1:
            aligner_LBC = LBC_aligner_dict[(len(L) // 10) * 10]
            LBC_match = aligner_LBC(L, clip_cigar=False)
            if (LBC_match.text_end - LBC_match.text_start)>= 6 and LBC_match.score > -14:             
                LBC_end = LBC_match.text_end
        else:
            LBC_match = [Match(findLBC, findLBC+len(LBC),0)]
            LBC_end = LBC_match[0].end    
        if not LBC_end:
            not_found += 1
            continue
        double = False
        findfull = L.find(fullRepeat)
        if findfull ==-1:
            aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
            drfull_match = aligner_drfull(L, clip_cigar=False)
            full_start = drfull_match.text_start
            full_end = drfull_match.text_end
            double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
        elif findfull !=-1:
            full_start = findfull
            full_end = findfull + len(fullRepeat)
            double = True

        if double:  # if full repeat is present
            # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
            tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
            if len(tempSpacerProximal) < 50:  
                L = L[len(tempSpacerProximal) - 5:]
                double = False
                if len(L) >= 60:
                    aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                    drfull_match = aligner_drfull(L, clip_cigar=False)
                    full_start = drfull_match.text_start
                    full_end = drfull_match.text_end
                    double = drfull_match.score > -16 and (drfull_match.text_end - drfull_match.text_start) >= 25
                    tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
                    if len(tempSpacerProximal) < 60:
                        double = False
                else:
                    continue
                
        if double:
            full_repeat += 1
            print("drfull match:")
            print(drfull_match.pretty)
            if  findLBC==-1:
                print("LBC match:")
                print(LBC_match.pretty)

    with open(txt_file, 'w') as file:
        sys.stdout = file
        print("tot: "+str(tot))
        print("not found: "+str(not_found))
        print("full repeats in LBC reads: "+str(full_repeat))
    sys.stdout = sys.__stdout__

    print("tot: "+str(tot))
    print("not found: "+str(not_found))
    print("full repeats in LBC reads: "+str(full_repeat))



################## EXPERIMENT 3: Here I test how many of the reads that have a LBC are double acquisitions (fullrepeat present) using the current fuzzysearch approach #########################
#
# Findings: the fuzzysearch allowed mismatches for both the full repeat and LBC are more strict than what I have with WFA2. Only 110 reads that have a LBC, also have a full repeat.
#           None of them are real Dr1-Spacer-Drfull-spacer-dr2-LBC constructs. They are weird reads where eg. the full repeat comes really early in the read.
#
if experiment_3:
    txt_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_3_{files[fasta_file]}"
    not_found = 0
    tot = 0
    full_repeat = 0
    for L in F:
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            continue
        
        L=L.strip()
        # ignore reads that are too short to detect adapted spacer
        if len(L) < 60: 
            continue
        tot += 1
        findLBC = L.find(LBC)
        if findLBC==-1:
            numMismatches = int(round(len(LBC)*0.1))
            LBCcoord = fuzzysearch.find_near_matches(LBC, L, max_l_dist = numMismatches)
        else:
            LBCcoord = [Match(findLBC, findLBC+len(LBC),0)]
        if not LBCcoord:
            not_found += 1
            continue
        numMismatches = 2
        tempFullRepeat = fuzzysearch.find_near_matches(fullRepeat, L, max_l_dist = numMismatches)
        if tempFullRepeat:  # if full repeat is found
            if len(tempFullRepeat) == 1:     
                full_repeat += 1
                print(L)
                print("drfull match:")
                print(tempFullRepeat[0].start,tempFullRepeat[0].end)
                print("LBC match:")
                print(LBCcoord[0].start,LBCcoord[0].end)
        
        
    with open(txt_file, 'w') as file:
        sys.stdout = file
        print("tot: "+str(tot))
        print("not found: "+str(not_found))
        print("full repeats in LBC reads: "+str(full_repeat))
    sys.stdout = sys.__stdout__
    print("tot: "+str(tot))
    print("not found: "+str(not_found))
    print("full repeats in LBC reads: "+str(full_repeat))


################## EXPERIMENT 4: Here I test if double acquisitions are unique a lot of the time (currently double ac. are not extracted because they lack the LBC as it is cut off) #########################
#
# Findings: 
# 
#
#
#
#

if experiment_4:
    txt_file = f"Experiment_4_{files[fasta_file]}"
    not_found = 0
    found = 0
    tot = 0
    none_extr = 0
    full_repeat = 0
    perfect_alignment = 0
    double_extr = 0
    proximal_extr = 0
    distal_extr = 0
    ds_already_there = 0
    close_ds_already_here = 0
    dists_from_ds_already_there = 0
    close_dists_from_ds_already_there = 0
    proxs_from_ds_already_there = 0
    close_proxs_from_ds_already_there = 0
    dists_from_SS_already_there = 0
    close_dists_from_SS_already_there = 0
    proxs_from_SS_already_there = 0
    close_proxs_from_SS_already_there = 0
    D = {}
    D_exact = {}
    D_1mm = {}
    D_2mm = {}
    D_3mm = {}
    D_4mm = {}
    counter = 0
    for L in F:
       counter += 1
       if counter > 100000:
            continue
       if counter % 10000 == 0:
           print(counter)
       if '>' in L: # defline, skip for processing but save read name
           readName = L.strip()
           continue
        
       L=L.strip()
       # ignore reads that are too short to detect adapted spacer
       if len(L) < 60: 
           continue
       tot += 1
    
       double = False
       findfull = L.find(fullRepeat)   
       if findfull ==-1:
           aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
           drfull_match = aligner_drfull(L, clip_cigar=False)
           full_start = drfull_match.text_start
           full_end = drfull_match.text_end
           double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
       elif findfull !=-1:
           full_start = findfull
           full_end = findfull + len(fullRepeat)
           double = True
           perfect_alignment += 1
    
       if double:  # if full repeat is present
           # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
           tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
           if len(tempSpacerProximal) < 60:  
               L = L[len(tempSpacerProximal) - 5:]
               double = False
               if len(L) >= 60:
                   aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                   drfull_match = aligner_drfull(L, clip_cigar=False)
                   full_start = drfull_match.text_start
                   full_end = drfull_match.text_end
                   double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
                   tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
                   if len(tempSpacerProximal) < 60:
                       double = False
               else:
                   continue
                
       if double:
               full_repeat += 1
               spacerProximal = ''
               spacerProximal = editSpacer_WFA(tempSpacerProximal,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict)
               spacerDistal = ''
               tempSpacerDistal = L[(full_end-len(firstRepeat)-8):] # the for is to account for the tendency of regex to add up to 3 nucleotides at the end of a spacer
               if len(tempSpacerDistal) >= 40:
                       spacerDistal = editSpacer_WFA(tempSpacerDistal,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict) 

               if spacerDistal and spacerProximal:
                   double_extr += 1
                   doubleSpacer= spacerProximal+"_"+spacerDistal                
                   if doubleSpacer in D:
                       ds_already_there += 1  
                       D[doubleSpacer][1]+=1
                   else:
                       aligner = WavefrontAligner(doubleSpacer, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                       present = 0
                       for key in D.keys():
                           match = aligner(key, clip_cigar=False)
                           if match.score > -30:
                               close_ds_already_here +=1
                               D[key][1]+=1
                               print(match.pretty)
                               print(match.score)
                               present = 1
                       if present == 0:
                               D[doubleSpacer]=[readName+'_doubleAcquisitions_both',0]
                               D[doubleSpacer][1]+=1
                   if spacerDistal in D:
                       dists_from_ds_already_there += 1
                       D[spacerDistal][1]+=1
                   else:
                       aligner = WavefrontAligner(spacerDistal, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                       present = 0
                       for key in D.keys():
                           match = aligner(key, clip_cigar=False)
                           if match.score > -30:
                               close_dists_from_ds_already_there +=1
                               D[key][1]+=1
                               present = 1
                       if present == 0:
                               D[spacerDistal]=[readName+'_doubleAcquisitions_both',0]
                               D[spacerDistal][1]+=1
                   if spacerProximal in D:
                       proxs_from_ds_already_there += 1
                       D[spacerProximal][1]+=1
                   else:
                       aligner = WavefrontAligner(spacerProximal, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                       present = 0
                       for key in D.keys():
                           match = aligner(key, clip_cigar=False)
                           if match.score > -30:
                               close_proxs_from_ds_already_there +=1
                               D[key][1]+=1
                               present = 1
                       if present == 0:
                               D[spacerProximal]=[readName+'_doubleAcquisitions_both',0]
                               D[spacerProximal][1]+=1
               elif spacerDistal:
                   distal_extr +=1
                   if spacerDistal in D:
                       D[spacerDistal][1]+=1
                       dists_from_SS_already_there += 1
                   else:
                       aligner = WavefrontAligner(spacerDistal, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                       present = 0
                       for key in D.keys():
                           match = aligner(key, clip_cigar=False)
                           if match.score > -30:
                               D[key][1]+=1
                               close_dists_from_SS_already_there +=1
                               present = 1
                       if present == 0:
                               D[spacerDistal]=[readName+'_doubleAcquisitions_both',0]
                               D[spacerDistal][1]+=1

               elif spacerProximal:
                   proximal_extr +=1
                   if spacerProximal in D:
                       D[spacerProximal][1]+=1
                       proxs_from_SS_already_there += 1
                   else:
                       aligner = WavefrontAligner(spacerProximal, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                       present = 0
                       for key in D.keys():
                           match = aligner(key, clip_cigar=False)
                           if match.score > -30:
                               close_proxs_from_SS_already_there +=1
                               D[key][1]+=1
                               present = 1
                       if present == 0:
                               D[spacerProximal]=[readName+'_doubleAcquisitions_both',0]
                               D[spacerProximal][1]+=1
               else:
                   none_extr += 1
       else:
           continue
    with open(txt_file, 'w') as file:
        sys.stdout = file
        print("tot reads: "+str(tot))
        print("full repeats totally: "+str(full_repeat))
        print("Full read found with no mismatches: "+str(perfect_alignment))
        print("full repeat found but no spacer could be extracted: "+str(none_extr))
        print("full repeat found and 2 spacers extracted: "+str(double_extr))
        print("full repeat found and only proximal spacer extracted: "+str(proximal_extr))
        print("full repeat found and only distal spacer extracted: "+str(distal_extr))
        print("double spacer extracted that has already been found prior without mismatches: "+str(ds_already_there))
        print("double spacer extracted that has already been found prior with some mismatches: "+str(close_ds_already_here))
        print("distal spacer extracted from double spacer that has already been found prior without mismatches: "+str(dists_from_ds_already_there))
        print("distal spacer extracted from double spacer that has already been found prior with some mismatches: "+str(close_dists_from_ds_already_there))
        print("proximal spacer extracted from double spacer that has already been found prior without mismatches: "+str(proxs_from_ds_already_there))
        print("proximal spacer extracted from double spacer that has already been found prior with some mismatches: "+str(close_proxs_from_ds_already_there))
        print("(only) distal spacer extracted from double ac. that has already been found prior without mismatches: "+str(dists_from_SS_already_there))
        print("(only) distal spacer extracted from double ac. that has already been found prior with some mismatches: "+str(close_dists_from_SS_already_there))
        print("(only) proximal spacer extracted from double ac. that has already been found prior without mismatches: "+str(proxs_from_SS_already_there))
        print("(only) proximal spacer extracted from double ac. that has already been found prior with some mismatches: "+str(close_proxs_from_SS_already_there))
    sys.stdout = sys.__stdout__    
    print("tot reads: "+str(tot))
    print("full repeats totally: "+str(full_repeat))
    print("Full read found with no mismatches: "+str(perfect_alignment))
    print("full repeat found but no spacer could be extracted: "+str(none_extr))
    print("full repeat found and 2 spacers extracted: "+str(double_extr))
    print("full repeat found and only proximal spacer extracted: "+str(proximal_extr))
    print("full repeat found and only distal spacer extracted: "+str(distal_extr))
    print("double spacer extracted that has already been found prior without mismatches: "+str(ds_already_there))
    print("double spacer extracted that has already been found prior with some mismatches: "+str(close_ds_already_here))
    print("distal spacer extracted from double spacer that has already been found prior without mismatches: "+str(dists_from_ds_already_there))
    print("distal spacer extracted from double spacer that has already been found prior with some mismatches: "+str(close_dists_from_ds_already_there))
    print("proximal spacer extracted from double spacer that has already been found prior without mismatches: "+str(proxs_from_ds_already_there))
    print("proximal spacer extracted from double spacer that has already been found prior with some mismatches: "+str(close_proxs_from_ds_already_there))
    print("(only) distal spacer extracted from double ac. that has already been found prior without mismatches: "+str(dists_from_SS_already_there))
    print("(only) distal spacer extracted from double ac. that has already been found prior with some mismatches: "+str(close_dists_from_SS_already_there))
    print("(only) proximal spacer extracted from double ac. that has already been found prior without mismatches: "+str(proxs_from_SS_already_there))
    print("(only) proximal spacer extracted from double ac. that has already been found prior with some mismatches: "+str(close_proxs_from_SS_already_there))




################## EXPERIMENT 5: Here I test if single acquisitions are unique a lot of the time  #########################
# I look at the first 1'000'000 reads in the fasta file
# Findings: 
# 
#
#
#
#
#
if experiment_5:
    txt_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_5_summary_2{files[fasta_file]}"
    full_repeat_found = 0
    tot = 0
    single_aq = 0
    not_identified = 0
    top_spacers = 1000
    D_exact = {}
    D_1mm = {}
    D_2mm = {}
    D_3mm = {}
    D_4mm = {}
    counter = 0
    for L in F:
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            continue
    
        L=L.strip()
        # ignore reads that are too short to detect adapted spacer
        counter += 1
        if counter > 10000:
             break
        if counter % 1000 == 0:
            print(counter)
        if len(L) < 60: 
            continue
        tot += 1
    
        double = False
        findfull = L.find(fullRepeat)   
        if findfull ==-1:
            aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
            drfull_match = aligner_drfull(L, clip_cigar=False)
            full_start = drfull_match.text_start
            full_end = drfull_match.text_end
            double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
        elif findfull !=-1:
            full_start = findfull
            full_end = findfull + len(fullRepeat)
            double = True
    
        if double:  # if full repeat is present
            # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
            tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
            if len(tempSpacerProximal) < 60:  
                L = L[len(tempSpacerProximal) - 5:]
                double = False
                if len(L) >= 60:
                    aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                    drfull_match = aligner_drfull(L, clip_cigar=False)
                    full_start = drfull_match.text_start
                    full_end = drfull_match.text_end
                    double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
                    tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
                    if len(tempSpacerProximal) < 60:
                        double = False
                else:
                    continue
    
        if double:
            full_repeat_found +=1
            continue
        else:
            single_aq += 1
            spacer = editSpacer_WFA(L,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict)    
            if not spacer:
                not_identified += 1
                continue      
            if spacer in D_exact:
                D_exact[spacer][1]+=1
                
            else:
                D_exact[spacer]=[readName,0]
                D_exact[spacer][1]+=1
            aligner = WavefrontAligner(spacer, gap_opening = 20, gap_opening2 = 24, gap_extension = 10, span="end-to-end")
            present = 0
            for key in D_1mm.keys():
                    match = aligner(key, clip_cigar=False)  #erst bester der matched. wenn es mehrere gibt, wird dies nicht festgestellt 
                    if match.score >= -4:
                        D_1mm[key][1]+=1
                        present = 1
                        break
            if present == 0:
                    D_1mm[spacer]=[readName,0]
                    D_1mm[spacer][1]+=1
            
            present = 0
            for key in D_2mm.keys():
                    match = aligner(key, clip_cigar=False)  #erst bester der matched. wenn es mehrere gibt, wird dies nicht festgestellt 
                    if match.score >= -8:
                        D_2mm[key][1]+=1
                        present = 1
                        break
            if present == 0:
                    D_2mm[spacer]=[readName,0]
                    D_2mm[spacer][1]+=1
            present = 0
            for key in D_3mm.keys():
                match = aligner(key, clip_cigar=False)  #erst bester der matched. wenn es mehrere gibt, wird dies nicht festgestellt 
                if match.score >= -12:
                    D_3mm[key][1]+=1
                    present = 1
                    break
            if present == 0:
                    D_3mm[spacer]=[readName,0]
                    D_3mm[spacer][1]+=1
            present = 0
            for key in D_4mm.keys():
                match = aligner(key, clip_cigar=False)  #erst bester der matched. wenn es mehrere gibt, wird dies nicht festgestellt 
                if match.score >= -16:
                    D_4mm[key][1]+=1
                    present = 1
                    break
            if present == 0:
                    D_4mm[spacer]=[readName,0]
                    D_4mm[spacer][1]+=1
            
    with open(txt_file, 'w') as file:
        sys.stdout = file
        print("tot reads: "+str(tot))
        print("Full repeat found: "+str(full_repeat_found))
        print("Single acquisition: "+str(single_aq))
        print("Single acquisition but no spacer identified: "+ str(not_identified))
        print("Unique spacers exact match:" + str(len(D_exact)))
        print("Unique spacers 1 mm allowed:" + str(len(D_1mm)))
        print("Unique spacers 2 mm allowed:" + str(len(D_2mm)))
        print("Unique spacers 3 mm allowed:" + str(len(D_3mm)))
        print("Unique spacers 4 mm allowed:" + str(len(D_4mm)))
    sys.stdout = sys.__stdout__
 
    print("tot reads: "+str(tot))
    print("Full repeat found: "+str(full_repeat_found))
    print("Single acquisition: "+str(single_aq))
    print("Single acquisition but no spacer identified: "+ str(not_identified))
    print("Unique spacers exact match:" + str(len(D_exact)))
    print("Unique spacers 1 mm allowed:" + str(len(D_1mm)))
    print("Unique spacers 2 mm allowed:" + str(len(D_2mm)))
    print("Unique spacers 3 mm allowed:" + str(len(D_3mm)))
    print("Unique spacers 4 mm allowed:" + str(len(D_4mm)))

    dictionaries = [D_exact,D_1mm,D_2mm,D_3mm,D_4mm]

    combined_data = []
    for index, dictionary in enumerate(dictionaries, start=1):
        sorted_items = sorted(dictionary.items(), key=lambda x: x[1][1], reverse=True)
        top_n_spacers = sorted_items[:top_spacers]
       

        for rank, (spacer, count) in enumerate(top_n_spacers, start=1):
            combined_data.append((f'D{index}', rank, spacer, count[1]))
    csv_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_5_2_{files[fasta_file]}.csv"
    with open(csv_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Dictionary', 'Rank', 'Spacer', 'Count'])
        csv_writer.writerows(combined_data)


################## EXPERIMENT 6: Full repeats how many uniques #########################
#
# Findings: 
# 
#
#
#
#
#
if experiment_6:
    txt_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_6_summary_{files[fasta_file]}"
    full_repeat = 0
    tot = 0
    proximal_already_seen_as_distal = 0
    proximal_already_seen_as_proximal = 0
    distal_already_seen_as_distal = 0
    distal_already_seen_as_proximal = 0
    perfect_alignment = 0
    double_extr = 0
    D_single_exact = {}
    D_double_exact = {}
    D_double_1mm = {}
    D_double_2mm = {}
    D_double_3mm = {}
    D_double_4mm = {}
    D_double_distal = {}
    D_double_proximal = {}
    D_proximal_in_single = {}
    D_distal_in_single = {}
    top_spacers = 100

    for L in F:
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            continue
   
        L=L.strip()
        # ignore reads that are too short to detect adapted spacer
        if len(L) < 60: 
            continue
        double = False
        findfull = L.find(fullRepeat)   
        if findfull ==-1:
            aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
            drfull_match = aligner_drfull(L, clip_cigar=False)
            full_start = drfull_match.text_start
            full_end = drfull_match.text_end
            double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
        elif findfull !=-1:
            full_start = findfull
            full_end = findfull + len(fullRepeat)
            double = True
        
        if double:  # if full repeat is present
            # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
            tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
            if len(tempSpacerProximal) < 60:  
                L = L[len(tempSpacerProximal) - 5:]
                double = False
                if len(L) >= 60:
                    aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                    drfull_match = aligner_drfull(L, clip_cigar=False)
                    full_start = drfull_match.text_start
                    full_end = drfull_match.text_end
                    double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
                    tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
                    if len(tempSpacerProximal) < 60:
                        double = False
                else:
                    continue
   
        if double:
            continue
        else:
            spacer = editSpacer_WFA(L,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict)    
            if not spacer:
                continue      
            if spacer in D_single_exact:
                D_single_exact[spacer][1]+=1
                
            else:
                D_single_exact[spacer]=[readName,0]
                D_single_exact[spacer][1]+=1 
       

    print("step one done ")
    
    F = open(fasta_directory+files[fasta_file], mode='r')
    counter = 0
    for L in F:
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            continue
   
        L=L.strip()
        # ignore reads that are too short to detect adapted spacer
        if len(L) < 60: 
            continue
        counter += 1
        if counter % 50000 == 0:
            print(counter)
        tot += 1
        double = False
        findfull = L.find(fullRepeat)   
        if findfull ==-1:
            aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
            drfull_match = aligner_drfull(L, clip_cigar=False)
            full_start = drfull_match.text_start
            full_end = drfull_match.text_end
            double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
        elif findfull !=-1:
            full_start = findfull
            full_end = findfull + len(fullRepeat)
            double = True
            perfect_alignment += 1
   
        if double:  # if full repeat is present
            # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
            tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
            if len(tempSpacerProximal) < 60:  
                L = L[len(tempSpacerProximal) - 5:]
                double = False
                if len(L) >= 60:
                    aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                    drfull_match = aligner_drfull(L, clip_cigar=False)
                    full_start = drfull_match.text_start
                    full_end = drfull_match.text_end
                    double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
                    tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
                    if len(tempSpacerProximal) < 60:
                        double = False
                else:
                    continue
   
        if double:
                full_repeat += 1
                spacerProximal = ''
                spacerProximal = editSpacer_WFA(tempSpacerProximal,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict)
                spacerDistal = ''
                tempSpacerDistal = L[(full_end-len(firstRepeat)-8):] # the for is to account for the tendency of regex to add up to 3 nucleotides at the end of a spacer
                if len(tempSpacerDistal) >= 40:
                        spacerDistal = editSpacer_WFA(tempSpacerDistal,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict) 
               
                if spacerDistal and spacerProximal:
                    double_extr += 1
                    doubleSpacer = spacerProximal+"_"+spacerDistal                
                    if doubleSpacer in D_double_exact:
                        D_double_exact[doubleSpacer][1]+=1
                
                    else:
                        D_double_exact[doubleSpacer]=[readName,0]
                        D_double_exact[doubleSpacer][1]+=1
                    aligner = WavefrontAligner(doubleSpacer, gap_opening = 20, gap_opening2 = 24, gap_extension = 10, span="end-to-end")
                    present = 0
                    for key in D_double_1mm.keys():
                            match = aligner(key, clip_cigar=False)  
                            if match.score >= -4:
                                D_double_1mm[key][1]+=1
                                present = 1
                                break
                    if present == 0:
                            D_double_1mm[doubleSpacer]=[readName,0]
                            D_double_1mm[doubleSpacer][1]+=1
            
                    present = 0
                    for key in D_double_2mm.keys():
                            match = aligner(key, clip_cigar=False)  
                            if match.score >= -8:
                                D_double_2mm[key][1]+=1
                                present = 1
                                break
                    if present == 0:
                            D_double_2mm[doubleSpacer]=[readName,0]
                            D_double_2mm[doubleSpacer][1]+=1
                    present = 0
                    for key in D_double_3mm.keys():
                        match = aligner(key, clip_cigar=False)   
                        if match.score >= -12:
                            D_double_3mm[key][1]+=1
                            present = 1
                            break
                    if present == 0:
                            D_double_3mm[doubleSpacer]=[readName,0]
                            D_double_3mm[doubleSpacer][1]+=1
                    present = 0
                    for key in D_double_4mm.keys():
                        match = aligner(key, clip_cigar=False)  
                        if match.score >= -16:
                            D_double_4mm[key][1]+=1
                            present = 1
                            break
                    if present == 0:
                            D_double_4mm[doubleSpacer]=[readName,0]
                            D_double_4mm[doubleSpacer][1]+=1
                    
                    if spacerDistal in D_single_exact:
                        if spacerDistal in D_distal_in_single:
                            D_distal_in_single[spacerDistal][1]+=1
                        else:
                            D_distal_in_single[spacerDistal]=[readName,0]
                            D_distal_in_single[spacerDistal][1]+=1
            
                    if spacerProximal in D_single_exact:
                        if spacerProximal in D_proximal_in_single:
                            D_proximal_in_single[spacerProximal][1]+=1
                        else:
                            D_proximal_in_single[spacerProximal]=[readName,0]
                            D_proximal_in_single[spacerProximal][1]+=1
                    
                    if spacerDistal in D_double_distal:
                        distal_already_seen_as_distal +=1
                        D_double_distal[spacerDistal][1]+=1
                    
                    else:
                        D_double_distal[spacerDistal]=[readName,0]
                        D_double_distal[spacerDistal][1]+=1
                    
                    if spacerProximal in D_double_proximal:
                        proximal_already_seen_as_proximal +=1
                        D_double_proximal[spacerProximal][1]+=1
                    
                    else:
                        D_double_proximal[spacerProximal]=[readName,0]
                        D_double_proximal[spacerProximal][1]+=1
                    
                    if spacerDistal in D_double_proximal:
                        distal_already_seen_as_proximal +=1
                    if spacerProximal in D_double_distal:
                        proximal_already_seen_as_distal +=1

        else:
            continue


    with open(txt_file, 'w') as file:
            sys.stdout = file
            print("tot reads: "+str(tot)) # yes
            print("full repeats totally: "+str(full_repeat)) # yes
            print("full repeats perfect alignment: "+str(perfect_alignment)) # yes
            print("full repeat found and 2 spacers extracted: "+str(double_extr)) # yes
            print("Distal spacer already seen as distal spacer in other double acquisition event: "+str(distal_already_seen_as_distal)) # yes
            print("Proximal spacer already seen as proximal spacer in other double acquisition event: "+str(proximal_already_seen_as_proximal)) # yes
            print("Distal spacer already seen as proximal spacer in other double acquisition event: "+str(distal_already_seen_as_proximal)) # yes
            print("Proximal spacer already seen as distal spacer in other double acquisition event: "+str(proximal_already_seen_as_distal)) # yes
            print("Number of unique double spacers exact match: " + str(len(D_double_exact)))
            print("Number of unique double spacers 1 mm allowed: " + str(len(D_double_1mm)))
            print("Number of unique double spacers 2 mm allowed: " + str(len(D_double_2mm)))
            print("Number of unique double spacers 3 mm allowed: " + str(len(D_double_3mm)))
            print("Number of unique double spacers 4 mm allowed: " + str(len(D_double_4mm)))
            print("Number of proximal spacers of double ac. events also identified as single ac.: " + str(len(D_proximal_in_single)))
            print("Number of distal spacers of double ac. events also identified as single ac.: " + str(len(D_distal_in_single)))
            print("Number of unique distal spacers (as distal spacers): "+str(len(D_double_distal))) # yes
            print("Number of unique proximal spacers (as proximal spacers): "+str(len(D_double_proximal))) # yes
    sys.stdout = sys.__stdout__

    print("tot reads: "+str(tot)) # yes
    print("full repeats totally: "+str(full_repeat)) # yes
    print("full repeats perfect alignment: "+str(perfect_alignment)) # yes
    print("full repeat found and 2 spacers extracted: "+str(double_extr)) # yes
    print("Distal spacer already seen as distal spacer in other double acquisition event: "+str(distal_already_seen_as_distal)) # yes
    print("Proximal spacer already seen as proximal spacer in other double acquisition event: "+str(proximal_already_seen_as_proximal)) # yes
    print("Distal spacer already seen as proximal spacer in other double acquisition event: "+str(distal_already_seen_as_proximal)) # yes
    print("Proximal spacer already seen as distal spacer in other double acquisition event: "+str(proximal_already_seen_as_distal)) # yes
    print("Number of unique double spacers exact match: " + str(len(D_double_exact)))
    print("Number of unique double spacers 1 mm allowed: " + str(len(D_double_1mm)))
    print("Number of unique double spacers 2 mm allowed: " + str(len(D_double_2mm)))
    print("Number of unique double spacers 3 mm allowed: " + str(len(D_double_3mm)))
    print("Number of unique double spacers 4 mm allowed: " + str(len(D_double_4mm)))
    print("Number of proximal spacers of double ac. events also identified as single ac.: " + str(len(D_proximal_in_single)))
    print("Number of distal spacers of double ac. events also identified as single ac.: " + str(len(D_distal_in_single)))
    print("Number of unique distal spacers (as distal spacers): "+str(len(D_double_distal))) # yes
    print("Number of unique proximal spacers (as proximal spacers): "+str(len(D_double_proximal))) # yes


    dictionaries = [D_single_exact,D_double_exact,D_double_1mm,D_double_2mm,D_double_3mm,D_double_4mm,D_double_distal,D_double_proximal,D_proximal_in_single,D_distal_in_single]
    dictionary_names = ["D_single_exact","D_double_exact","D_double_1mm","D_double_2mm","D_double_3mm","D_double_4mm","D_double_distal","D_double_proximal","D_proximal_in_single","D_distal_in_single"]
    combined_data = []
    for index, dictionary in enumerate(dictionaries, start=1):
        sorted_items = sorted(dictionary.items(), key=lambda x: x[1][1], reverse=True)
        top_n_spacers = sorted_items[:top_spacers]

        for rank, (spacer, count) in enumerate(top_n_spacers, start=1):
            combined_data.append((dictionary_names[index-1], rank, spacer, count[1])) #f'D{index}'
    csv_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_6_{files[fasta_file]}.csv"
    with open(csv_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Dictionary', 'Rank', 'Spacer', 'Count'])
        csv_writer.writerows(combined_data)

################## EXPERIMENT 6.5: Full repeats how many uniques #########################
#
# Findings: 
# 
#
#
#
#
#
if experiment_6_5:
    txt_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_6.5_summary_{files[fasta_file]}"
    full_repeat = 0
    tot = 0
    proximal_already_seen_as_distal = 0
    proximal_already_seen_as_proximal = 0
    distal_already_seen_as_distal = 0
    distal_already_seen_as_proximal = 0
    perfect_alignment = 0
    double_extr = 0
    prox_extr = 0
    D_single_exact = {}
    D_double_singles_exact = {}
    D_double_singles_1mm = {}
    D_double_distal = {}
    D_double_proximal = {}
    D_proximal_in_single = {}
    D_distal_in_single = {}
    top_spacers = 100

    for L in F:
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            continue
   
        L=L.strip()
        # ignore reads that are too short to detect adapted spacer
        if len(L) < 60: 
            continue
        double = False
        findfull = L.find(fullRepeat)   
        if findfull ==-1:
            aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
            drfull_match = aligner_drfull(L, clip_cigar=False)
            full_start = drfull_match.text_start
            full_end = drfull_match.text_end
            double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
        elif findfull !=-1:
            full_start = findfull
            full_end = findfull + len(fullRepeat)
            double = True
        
        if double:  # if full repeat is present
            # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
            tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
            if len(tempSpacerProximal) < 60:  
                L = L[len(tempSpacerProximal) - 5:]
                double = False
                if len(L) >= 60:
                    aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                    drfull_match = aligner_drfull(L, clip_cigar=False)
                    full_start = drfull_match.text_start
                    full_end = drfull_match.text_end
                    double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
                    tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
                    if len(tempSpacerProximal) < 60:
                        double = False
                else:
                    continue
   
        if double:
            continue
        else:
            spacer = editSpacer_WFA(L,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict)    
            if not spacer:
                continue      
            if spacer in D_single_exact:
                D_single_exact[spacer][1]+=1
                
            else:
                D_single_exact[spacer]=[readName,0]
                D_single_exact[spacer][1]+=1 
       

    print("step one done ")
    
    F = open(fasta_directory+files[fasta_file], mode='r')
    counter = 0
    for L in F:
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            continue
   
        L=L.strip()
        # ignore reads that are too short to detect adapted spacer
        if len(L) < 60: 
            continue
        counter += 1
        if counter % 50000 == 0:
            print(counter)
        tot += 1
        double = False
        findfull = L.find(fullRepeat)   
        if findfull ==-1:
            aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
            drfull_match = aligner_drfull(L, clip_cigar=False)
            full_start = drfull_match.text_start
            full_end = drfull_match.text_end
            double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
        elif findfull !=-1:
            full_start = findfull
            full_end = findfull + len(fullRepeat)
            double = True
            perfect_alignment += 1
   
        if double:  # if full repeat is present
            # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
            tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
            if len(tempSpacerProximal) < 60:  
                L = L[len(tempSpacerProximal) - 5:]
                double = False
                if len(L) >= 60:
                    aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                    drfull_match = aligner_drfull(L, clip_cigar=False)
                    full_start = drfull_match.text_start
                    full_end = drfull_match.text_end
                    double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
                    tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
                    if len(tempSpacerProximal) < 60:
                        double = False
                else:
                    continue
   
        if double:
                full_repeat += 1
                spacerProximal = ''
                spacerProximal = editSpacer_WFA(tempSpacerProximal,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict)
                spacerDistal = ''
                tempSpacerDistal = L[(full_end-len(firstRepeat)-8):] # the for is to account for the tendency of regex to add up to 3 nucleotides at the end of a spacer
                if len(tempSpacerDistal) >= 40:
                        spacerDistal = editSpacer_WFA(tempSpacerDistal,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict) 
               
                if spacerDistal and spacerProximal:
                    double_extr += 1 
                    aligner_proximal = WavefrontAligner(spacerProximal, gap_opening = 20, gap_opening2 = 24, gap_extension = 10, span="end-to-end")
                    aligner_distal = WavefrontAligner(spacerDistal, gap_opening = 20, gap_opening2 = 24, gap_extension = 10, span="end-to-end")
                    
                    if spacerDistal in D_double_singles_exact: 
                        D_double_singles_exact[spacerDistal][1]+=1
                    else:
                        D_double_singles_exact[spacerDistal]=[readName,0]
                        D_double_singles_exact[spacerDistal][1]+=1
                    
                    present = 0
                    for key in D_double_singles_1mm.keys():
                        match = aligner_distal(key, clip_cigar=False)  
                        if match.score >= -4:
                            D_double_singles_1mm[key][1]+=1
                            present = 1
                            break
                    if present == 0:
                            D_double_singles_1mm[spacerDistal]=[readName,0]
                            D_double_singles_1mm[spacerDistal][1]+=1

                    if spacerProximal in D_double_singles_exact: 
                        D_double_singles_exact[spacerProximal][1]+=1
                    else:
                        D_double_singles_exact[spacerProximal]=[readName,0]
                        D_double_singles_exact[spacerProximal][1]+=1  

                    present = 0
                    for key in D_double_singles_1mm.keys():
                        match = aligner_proximal(key, clip_cigar=False)  
                        if match.score >= -4:
                            D_double_singles_1mm[key][1]+=1
                            present = 1
                            break
                    if present == 0:
                            D_double_singles_1mm[spacerProximal]=[readName,0]
                            D_double_singles_1mm[spacerProximal][1]+=1

                    if spacerDistal in D_single_exact: 
                        if spacerDistal in D_distal_in_single:
                            D_distal_in_single[spacerDistal][1]+=1
                        else:
                            D_distal_in_single[spacerDistal]=[readName,0]
                            D_distal_in_single[spacerDistal][1]+=1
            
                    if spacerProximal in D_single_exact:
                        if spacerProximal in D_proximal_in_single:
                            D_proximal_in_single[spacerProximal][1]+=1
                        else:
                            D_proximal_in_single[spacerProximal]=[readName,0]
                            D_proximal_in_single[spacerProximal][1]+=1
                    
                    if spacerDistal in D_double_distal:
                        distal_already_seen_as_distal +=1
                        D_double_distal[spacerDistal][1]+=1
                    
                    else:
                        D_double_distal[spacerDistal]=[readName,0]
                        D_double_distal[spacerDistal][1]+=1
                    
                    if spacerProximal in D_double_proximal:
                        proximal_already_seen_as_proximal +=1
                        D_double_proximal[spacerProximal][1]+=1
                    
                    else:
                        D_double_proximal[spacerProximal]=[readName,0]
                        D_double_proximal[spacerProximal][1]+=1
                    
                    if spacerDistal in D_double_proximal:
                        distal_already_seen_as_proximal +=1
                    if spacerProximal in D_double_distal:
                        proximal_already_seen_as_distal +=1
                
                elif spacerProximal:
                    prox_extr += 1
                    aligner_proximal = WavefrontAligner(spacerProximal, gap_opening = 20, gap_opening2 = 24, gap_extension = 10, span="end-to-end")
                    if spacerProximal in D_double_singles_exact: 
                        D_double_singles_exact[spacerProximal][1]+=1
                    else:
                        D_double_singles_exact[spacerProximal]=[readName,0]
                        D_double_singles_exact[spacerProximal][1]+=1

                    for key in D_double_singles_1mm.keys():
                        match = aligner_proximal(key, clip_cigar=False)  
                        if match.score >= -4:
                            D_double_singles_1mm[key][1]+=1
                            present = 1
                            break
                    if present == 0:
                            D_double_singles_1mm[spacerProximal]=[readName,0]
                            D_double_singles_1mm[spacerProximal][1]+=1

        else:
            continue


    with open(txt_file, 'w') as file:
            sys.stdout = file
            print("tot reads: "+str(tot)) # yes
            print("full repeats totally: "+str(full_repeat)) # yes
            print("full repeats perfect alignment: "+str(perfect_alignment)) # yes
            print("full repeat found and 2 spacers extracted: "+str(double_extr)) # yes
            print("full repeat found and only proximal spacer extracted: "+str(prox_extr))
            print("Number of spacers extracted from double acquisition events: " + str((double_extr * 2 + prox_extr)))
            print("Number of unique spacers exact match: " + str(len(D_double_singles_exact)))
            print("Number of unique spacers 1 mm allowed: " + str(len(D_double_singles_1mm)))
            print("Distal spacer already seen as distal spacer in other double acquisition event: "+str(distal_already_seen_as_distal)) # yes
            print("Proximal spacer already seen as proximal spacer in other double acquisition event: "+str(proximal_already_seen_as_proximal)) # yes
            print("Distal spacer already seen as proximal spacer in other double acquisition event: "+str(distal_already_seen_as_proximal)) # yes
            print("Proximal spacer already seen as distal spacer in other double acquisition event: "+str(proximal_already_seen_as_distal)) # yes
            print("Number of proximal spacers of double ac. events also identified as single ac.: " + str(len(D_proximal_in_single)))
            print("Number of distal spacers of double ac. events also identified as single ac.: " + str(len(D_distal_in_single)))
            print("Number of unique distal spacers (as distal spacers): "+str(len(D_double_distal))) # yes
            print("Number of unique proximal spacers (as proximal spacers): "+str(len(D_double_proximal))) # yes
    sys.stdout = sys.__stdout__

    print("tot reads: "+str(tot)) # yes
    print("full repeats totally: "+str(full_repeat)) # yes
    print("full repeats perfect alignment: "+str(perfect_alignment)) # yes
    print("full repeat found and 2 spacers extracted: "+str(double_extr)) # yes
    print("full repeat found and only proximal spacer extracted: "+str(prox_extr))
    print("Number of spacers extracted from double acquisition events: " + str((double_extr * 2 + prox_extr)))
    print("Number of unique spacers exact match: " + str(len(D_double_singles_exact)))
    print("Number of unique spacers 1 mm allowed: " + str(len(D_double_singles_1mm)))
    print("Distal spacer already seen as distal spacer in other double acquisition event: "+str(distal_already_seen_as_distal)) # yes
    print("Proximal spacer already seen as proximal spacer in other double acquisition event: "+str(proximal_already_seen_as_proximal)) # yes
    print("Distal spacer already seen as proximal spacer in other double acquisition event: "+str(distal_already_seen_as_proximal)) # yes
    print("Proximal spacer already seen as distal spacer in other double acquisition event: "+str(proximal_already_seen_as_distal)) # yes
    print("Number of proximal spacers of double ac. events also identified as single ac.: " + str(len(D_proximal_in_single)))
    print("Number of distal spacers of double ac. events also identified as single ac.: " + str(len(D_distal_in_single)))
    print("Number of unique distal spacers (as distal spacers): "+str(len(D_double_distal))) # yes
    print("Number of unique proximal spacers (as proximal spacers): "+str(len(D_double_proximal))) # yes

    dictionaries = [D_single_exact,D_double_singles_exact,D_double_singles_1mm,D_double_distal,D_double_proximal,D_proximal_in_single,D_distal_in_single]
    dictionary_names = ["D_single_exact","D_double_singles_exact","D_double_singles_1mm","D_double_distal","D_double_proximal","D_proximal_in_single","D_distal_in_single"]
    combined_data = []
    for index, dictionary in enumerate(dictionaries, start=1):
        sorted_items = sorted(dictionary.items(), key=lambda x: x[1][1], reverse=True)
        top_n_spacers = sorted_items[:top_spacers]

        for rank, (spacer, count) in enumerate(top_n_spacers, start=1):
            combined_data.append((dictionary_names[index-1], rank, spacer, count[1])) #f'D{index}'
    csv_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_6.5_{files[fasta_file]}.csv"
    with open(csv_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Dictionary', 'Rank', 'Spacer', 'Count'])
        csv_writer.writerows(combined_data)

################## EXPERIMENT 7: Analysing how much of the LBC remains at the right end of the read in the case of double acquisitions #########################
#
# Findings: 
# 
#
#
#
#
#
if experiment_7:
    txt_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_7_summary_{files[fasta_file]}"
    
    mean = 0
    M2 = 0 
    count = 0
    counter = 0
    double_identified = 0
    double_extracted = 0
    LBC_residuals = 0
    residual_length = []
    residual_match = []

    for L in F:
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            continue
   
        L=L.strip()
        # ignore reads that are too short to detect adapted spacer
        if len(L) < 60: 
            continue
        counter += 1
        if counter % 50000 == 0:
            print(counter)
        double = False
        findfull = L.find(fullRepeat)   
        if findfull ==-1:
            aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
            drfull_match = aligner_drfull(L, clip_cigar=False)
            full_start = drfull_match.text_start
            full_end = drfull_match.text_end
            double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
        elif findfull !=-1:
            full_start = findfull
            full_end = findfull + len(fullRepeat)
            double = True
   
        if double:  # if full repeat is present
            # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
            tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
            if len(tempSpacerProximal) < 60:  
                L = L[len(tempSpacerProximal) - 5:]
                double = False
                if len(L) >= 60:
                    aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                    drfull_match = aligner_drfull(L, clip_cigar=False)
                    full_start = drfull_match.text_start
                    full_end = drfull_match.text_end
                    double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
                    tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
                    if len(tempSpacerProximal) < 60:
                        double = False
                else:
                    continue
   
        if double:
                double_identified += 1
                spacerProximal = ''
                spacerProximal = editSpacer_WFA(tempSpacerProximal,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict)
                spacerDistal = ''
                tempSpacerDistal = L[(full_end-len(firstRepeat)-8):] # the for is to account for the tendency of regex to add up to 3 nucleotides at the end of a spacer
                if len(tempSpacerDistal) >= 40:
                        spacerDistal = editSpacer_WFA(tempSpacerDistal,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict) 
               
                if spacerDistal and spacerProximal:
                    double_extracted += 1
                    aligner = dr2_aligner_dict[(len(tempSpacerDistal) // 10) * 10]
                    match = aligner(tempSpacerDistal, clip_cigar=False)
                    x = len(tempSpacerDistal[match.text_end:])
                    x = x-9
                    if x > 0:
                        LBC_residuals += 1
                        residual_length.append(x)
                        if tempSpacerDistal[(match.text_end+9):(match.text_end+x+9)] == LBC[:x]:
                            residual_match.append(1)
                        else:
                            residual_match.append(0)
                    
                    count += 1
                    delta = x - mean
                    mean += delta / count
                    delta2 = x - mean
                    M2 += delta * delta2
        else:
            continue
    
    final_mean = mean
    variance_n = M2 / count
    final_stddev = math.sqrt(variance_n)

    with open(txt_file, 'w') as file:
        sys.stdout = file
        print("Reads processed: "+str(counter))
        print("Full repeats found: " + str(double_identified))
        print("Two spacers extracted: " + str(double_extracted))
        print("LBC residual present: " + str(LBC_residuals))
        print("LBC residual matches LBC: " + str(sum(residual_match)))
        print("Mean length of LBC residual: "+str(final_mean))
        print("Standard deviation of LBC residual length: "+str(final_stddev))

    sys.stdout = sys.__stdout__
    print("Reads processed: "+str(counter))
    print("Full repeats found: " + str(double_identified))
    print("Two spacers extracted: " + str(double_extracted))
    print("LBC residual present: " + str(LBC_residuals))
    print("LBC residual matches LBC: " + str(sum(residual_match)))
    print("Mean length of LBC residual: "+str(final_mean))
    print("Standard deviation of LBC residual length: "+str(final_stddev))
    
    data = list(zip(residual_length, residual_match))
    csv_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_7_LBC_residuals_{files[fasta_file]}.csv"
    with open(csv_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Residual_length', 'match'])
        csv_writer.writerows(data)

    
    
    
################## EXPERIMENT 8: Single acquisitions how many uniques #########################
#
# Findings: 
# 
#
#
#
#
#
if experiment_8:
    file_counter = 0
    for fasta_file in fasta_files:
        print("file: "+ str(file_counter))
        file_counter += 1
        txt_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_8_summary_{files[fasta_file]}"
        F = open(fasta_directory+files[fasta_file], mode='r')
        tot = 0
        single_aq = 0
        single_extract = 0
        D_single_exact = {}
        top_spacers = 1000
        for L in F:
            if '>' in L: # defline, skip for processing but save read name
                readName = L.strip()
                continue
            
            L=L.strip()
            # ignore reads that are too short to detect adapted spacer
            if len(L) < 60: 
                continue
            tot += 1
            if tot % 50000 == 0:
                print(tot)
            double = False
            findfull = L.find(fullRepeat)   
            if findfull ==-1:
                aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                drfull_match = aligner_drfull(L, clip_cigar=False)
                full_start = drfull_match.text_start
                full_end = drfull_match.text_end
                double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
            elif findfull !=-1:
                full_start = findfull
                full_end = findfull + len(fullRepeat)
                double = True

            if double:  # if full repeat is present
                # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
                tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
                if len(tempSpacerProximal) < 60:  
                    L = L[len(tempSpacerProximal) - 5:]
                    double = False
                    if len(L) >= 60:
                        aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                        drfull_match = aligner_drfull(L, clip_cigar=False)
                        full_start = drfull_match.text_start
                        full_end = drfull_match.text_end
                        double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
                        tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
                        if len(tempSpacerProximal) < 60:
                            double = False
                    else:
                        continue
                    
            if double:
                continue
            else:
                single_aq += 1
                spacer = editSpacer_WFA(L,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict)    
                if not spacer:
                    continue      
                single_extract += 1
                if spacer in D_single_exact:
                    D_single_exact[spacer][1]+=1

                else:
                    D_single_exact[spacer]=[readName,0]
                    D_single_exact[spacer][1]+=1 

        with open(txt_file, 'w') as file:
            sys.stdout = file
            print("tot reads: "+str(tot)) 
            print("No double aq: "+str(single_aq)) 
            print("Single spacer extracted: "+str(single_extract)) 
            print("Number of unique spacers: " + str(len(D_single_exact)))
            
        sys.stdout = sys.__stdout__
        dictionaries = [D_single_exact]
        dictionary_names = ["D_single_exact"]
        combined_data = []
        for index, dictionary in enumerate(dictionaries, start=1):
            sorted_items = sorted(dictionary.items(), key=lambda x: x[1][1], reverse=True)
            top_n_spacers = sorted_items[:top_spacers]
    
            for rank, (spacer, count) in enumerate(top_n_spacers, start=1):
                combined_data.append((dictionary_names[index-1], rank, spacer, count[1])) 
        csv_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_8_{files[fasta_file]}.csv"
        with open(csv_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(['Dictionary', 'Rank', 'Spacer', 'Count'])
            csv_writer.writerows(combined_data)

################## EXPERIMENT 8.5: Single acquisitions how many uniques less reads processed #########################
#
# Findings: 
# 
#
#
#
#
#
if experiment_8_5:
    file_counter = 0
    for fasta_file in fasta_files:
        print("file: "+ str(file_counter))
        file_counter += 1
        txt_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_8.5_summary_{files[fasta_file]}"
        F = open(fasta_directory+files[fasta_file], mode='r')
        tot = 0
        single_aq = 0
        single_extract = 0
        D_single_exact = {}
        top_spacers = 1000
        for L in F:
            if '>' in L: # defline, skip for processing but save read name
                readName = L.strip()
                continue
            
            L=L.strip()
            # ignore reads that are too short to detect adapted spacer
            if len(L) < 60: 
                continue
            tot += 1
            if tot % 50000 == 0:
                print(tot)
            if tot > 300000:
                break
            double = False
            findfull = L.find(fullRepeat)   
            if findfull ==-1:
                aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                drfull_match = aligner_drfull(L, clip_cigar=False)
                full_start = drfull_match.text_start
                full_end = drfull_match.text_end
                double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
            elif findfull !=-1:
                full_start = findfull
                full_end = findfull + len(fullRepeat)
                double = True

            if double:  # if full repeat is present
                # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
                tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
                if len(tempSpacerProximal) < 60:  
                    L = L[len(tempSpacerProximal) - 5:]
                    double = False
                    if len(L) >= 60:
                        aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                        drfull_match = aligner_drfull(L, clip_cigar=False)
                        full_start = drfull_match.text_start
                        full_end = drfull_match.text_end
                        double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
                        tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
                        if len(tempSpacerProximal) < 60:
                            double = False
                    else:
                        continue
                    
            if double:
                continue
            else:
                single_aq += 1
                spacer = editSpacer_WFA(L,firstRepeat,secondRepeat,20,60, dr1_aligner_dict, dr2_aligner_dict)    
                if not spacer:
                    continue      
                single_extract += 1
                if spacer in D_single_exact:
                    D_single_exact[spacer][1]+=1

                else:
                    D_single_exact[spacer]=[readName,0]
                    D_single_exact[spacer][1]+=1 

        with open(txt_file, 'w') as file:
            sys.stdout = file
            print("tot reads: "+str(tot)) 
            print("No double aq: "+str(single_aq)) 
            print("Single spacer extracted: "+str(single_extract)) 
            print("Number of unique spacers: " + str(len(D_single_exact)))
            
        sys.stdout = sys.__stdout__
        dictionaries = [D_single_exact]
        dictionary_names = ["D_single_exact"]
        combined_data = []
        for index, dictionary in enumerate(dictionaries, start=1):
            sorted_items = sorted(dictionary.items(), key=lambda x: x[1][1], reverse=True)
            top_n_spacers = sorted_items[:top_spacers]
    
            for rank, (spacer, count) in enumerate(top_n_spacers, start=1):
                combined_data.append((dictionary_names[index-1], rank, spacer, count[1])) 
        csv_file = f"/cluster/home/hugifl/recordseq-workflow-dev/dev-hugi/read_analysis_outputs/Experiment_8.5_{files[fasta_file]}.csv"
        with open(csv_file, 'w', newline='') as csvfile:
            csv_writer = csv.writer(csvfile)
            csv_writer.writerow(['Dictionary', 'Rank', 'Spacer', 'Count'])
            csv_writer.writerows(combined_data)