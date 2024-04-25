# Oct 13, 2023
# Florian Hugi, updated script from Tanmay Tanna and Wiona Glaenzer
# Script to extract spacers from reads.
# Changes from last version:
# - Wavefront alignment 2 (WFA2) is used to align the direct repeats to the reads --> faster, better performance (less FN)
# - In case of double acquisitions, the UMI and LBC are not required since they are cut off anyways (reads to short) --> spacers from double acquisition events are now extracted


from __future__ import division
import sys, os, argparse, operator, fuzzysearch
from collections import Counter
import collections
from editSpacer_WFA import editSpacer_WFA
from WFA_aligner_construction import construct_aligners
from pywfa import WavefrontAligner
from check_double import check_double

## set up parser for user inputs 
parser = argparse.ArgumentParser()

optional = parser._action_groups.pop()
required = parser.add_argument_group('required arguments')

## user inputs required
required.add_argument('-p', '--inFile', help='path to .fasta file.', dest='file_inFile')
required.add_argument('-o', '--outPath', help='path to output directory.', dest='path_outPath')
required.add_argument('-d1', '--drOne', help='first partial CRISPR repeat ', dest='dr_sequence1')
required.add_argument('-d2', '--drTwo', help='second partial CRISPR repeat ', dest='dr_sequence2')
required.add_argument('-df', '--drFull', help='second partial CRISPR repeat ', dest='dr_sequenceFull')



## optional (defaults provided)
optional.add_argument('-u', '--UMI', help='UMIs set on or off', dest='UMI_status', default='off')
optional.add_argument('-l', '--LBC', help='library identifier ', dest='LBC', default='null')
optional.add_argument('-t', '--UMIlength', help='length of UMIs', dest='UMI_length', type=int, default=10)
optional.add_argument('-s', '--outName', help='path to output directory.', dest='file_outName', default='null')
optional.add_argument('-a', help='Number of allowed mismatches in the first partial CRISPR repeat. Default=2', type=int, dest='m1', default=2)
optional.add_argument('-b', help='Number of allowed mismatches in the second partial CRISPR repeat. Default=0', type=int, dest='m2', default=0)
optional.add_argument('-c', help='Number of allowed mismatches in a single spacer sequence for collapsing if UMI sequence is identical. Default=4', type=int, dest='um1', default=4)
optional.add_argument('-d', help='Number of allowed mismatches in a multiple spacer sequence for collapsing if UMI sequence is identical. Default=6', type=int, dest='um2', default=6)
optional.add_argument('-m', help='Minimum spacer size. Default=20', type=int, dest='min', default=20)
optional.add_argument('-n', help='Maximum spacer size. Default=60', type=int, dest='max', default=60)
optional.add_argument('-sMin', help='Minimum stagger length according to primer design. Default=0', type=int, dest='min_stagger', default=0)
optional.add_argument('-sMax', help='Maximum stagger length according to primer design. Default=8', type=int, dest='max_stagger', default=8)
optional.add_argument('--infoFile', help='boolean to generate files with GC content and length info for spacers', dest='infoFile', action='store_true')
optional.add_argument('--no-infoFile', help='boolean to generate files with GC content and length info for spacers', dest='infoFile', action='store_false')
optional.add_argument('-co', '--contaminationstring', help='string with contamination sequences to quantify', default=None)                           ###### contaminations
optional.set_defaults(infoFile=True)

parser._action_groups.append(optional) 
args = parser.parse_args()

# assign arguments to variables

DR1Mismatch = int(args.m1)
DR2Mismatch = int(args.m2)
singleSpacerMismatch = int(args.um1)
multipleSpacerMismatch = int(args.um2)
minSpacer = int(args.min)
maxSpacer = int(args.max)
inFile = str(args.file_inFile)
outPath = str(args.path_outPath)+'/'
outName = str(args.file_outName)
firstRepeat = str(args.dr_sequence1)
secondRepeat = str(args.dr_sequence2)
fullRepeat = str(args.dr_sequenceFull)
UMI_status = str(args.UMI_status)
LBC = str(args.LBC)
UMI_length = int(args.UMI_length)
minStagger = int(args.min_stagger)
maxStagger = int(args.max_stagger)
infoFile=args.infoFile
if outName == 'null':
    outName = inFile.split("/")[-1]
if outName.endswith('.fasta'):
    outName = outName[:-6]
contaminations = str(args.contaminationstring)                                        ###### contaminations 
if contaminations is not None:
    if "," in contaminations:                                                           
        contaminationlist = contaminations.split(",")
        contaminationlist = list(map(str.strip, contaminationlist))
    else:
        contaminationlist = [contaminations]
    contamination_dict = {key: 0 for key in contaminationlist}


# define a custom data type for found matches
Match = collections.namedtuple('Match', ['start', 'end', 'dist'])

# construct the wfa2 alignment objects 
dr1_aligner_dict = construct_aligners(firstRepeat, 20, 150)  
dr2_aligner_dict = construct_aligners(secondRepeat, 20, 150) 
drfull_aligner_dict = construct_aligners(fullRepeat, 20, 150)
LBC_aligner_dict = construct_aligners(LBC, 20, 150)
 # identify and process files with the terms below
if ('.fasta' in inFile):
    
    # open inFile for reading/writing and report file being processed
    F = open(inFile, mode='r')
    G = open(outPath+outName+'.unique.fasta',mode='w')  # unique spacers based on spacer sequence only
    I = open(outPath+outName+'.doubleAcquisitions.fasta',mode='w') 
    # J = open(outPath+outName+'.doubleAcquisitions.paired.fasta',mode='w') # double acquisitions with both spacers 
    K = open(outPath+outName+'.all.fasta',mode='w')
    # MC = open(outPath+outName+'.multipleAcquisitions.complete.fasta',mode='w') # multiple acquisitions with all spacers
    # M = open(outPath+outName+'.multipleAcquisitions.fasta',mode='w')
    
    SS = open(outPath+"summaryStats.txt", mode = 'a')                 ##### tanmay commmit
    if contaminations is not None:
        CS = open(outPath+"contaminationStats.txt", mode = 'a')                                                                           ###### contaminations
    if infoFile:
        GC = open(outPath+outName+'.info.txt',mode='w')
        GC.write("Unique_Spacer_Sequence"+'\t'+"Sequence_Length"+'\t'+"GC_content"+'\t'+"count"+'\n')

    if UMI_status =='on':
        U = open(outPath+outName+'.umi.fasta',mode='w') # unique spacers based on unique combinations of spacer and umi sequence
        # UR = open(outPath+outName+'.umilist.txt',mode='w') # list of unique umis with the number of unique spacers associated to them 
        US = open(outPath+"umiStats.txt", mode = 'a')                                                                   ##### tanmay commit
        # UR.write("UMI"  + '\t' "unique_spacer_sequence" + '\t' "amplicon_number"+ '\n')
        if infoFile:
            UGC = open(outPath+outName+'.umi.info.txt',mode='w')
            UGC.write("Unique_Spacer_Sequence"+'\t'+"UMI"+'\t'+"Sequence_Length"+'\t'+"GC_content"+'\t'+"count"+'\n')



    os.system(str("echo '##################################################'"))
    os.system(str('echo '+"'"+inFile+' accepted for processing'+"'"))
    os.system(str('echo '+"'"+' UMI based unique spacers : '+UMI_status+"'"))

    readName = ''
    D={}
    rawReads=0 # total reads in fasta
    spacerReads=0 # reads having a spacer
    UniqueSingleAcquisitions=0 # number of unique single acquisitions based on spacer sequence
    UniqueDoubleAcquisitions=0 # number of unique double acquisitions based on spacer sequence
    # UniqueMultipleAcquisitions=0 # number of unique multiple acquisitions based on spacer sequence
    SinglefullRepeatReads=0 # reads with a double acquisition
    # MultiplefullRepeatReads=0 # reads with multiple acquisitions
    spacerReadsDoubleBoth=0 # double acquisitions with both spacers
    spacerReadsDoubleOne=0 # double acquisitions with one spacer
    # spacerReadsMultiComplete=0 # multiple acquisition with all spacers
    # spacerReadsMultiSome=0 # multiple acquisiton with one spacer
    # spacerReadsMultiNoSpacerBetweenFullDRs=0 # multiple DRs without spacers between them (probably PCR artifacts)
    minReadLength=minStagger + len(firstRepeat) + minSpacer + len(secondRepeat)
    allDistalSpacers=[]
    if 'null' not in LBC:
        NonLBCReads=0

    ## dictionaries and counters for UMI 

    if UMI_status =='on':
        UMISpacers={} # dictionary of UMIs with corresponding spacers as well as number of amplicons detected for spacer-umi combination
        # UMIcounts={} # dictionary of UMIs with number of associated unique spacers
        # products having only distal spacers of double/multiple acquisitions by using the UMI.
        UniqueSpacersUMI={} # dictionary of unique spacers with associated UMIs
        UniqueSpacersUMIReadNames={} # dictionary of read names corresponding to the above dictionary
        UniqueAcquisitionsUMI=0
        UniqueUMIs=0 # number of unique UMIs
        NonUMIReads=0 # number of reads without UMI upstream sequence
        UniqueSingleAcquisitionsUMI=0 # number of unique single acquisitions based on spacer and umi sequence
        UniqueDoubleAcquisitionsUMI=0 # number of unique double acquisitions based on spacer and umi sequence
        # UniqueMultipleAcquisitionsUMI=0 # number of unique multiple acquisitions based on spacer and umi sequence

    for L in F: # loop through reads in file
        if '>' in L: # defline, skip for processing but save read name
            readName = L.strip()
            rawReads+=1
            continue
        L=L.strip()
            
        # ignore reads that are too short to detect adapted spacer
        if len(L) < minReadLength: 
            continue
        
        skip_read, double, tempSpacerProximal, full_end, L = check_double(L, fullRepeat, drfull_aligner_dict, minReadLength, secondRepeat)
            
        if skip_read:
            if contaminations is not None:                                                                                                                            ###### contaminations
                    for contamination in contaminationlist:
                        if contamination in L:
                            contamination_dict[contamination] += 1
            continue
            

        # Identify and extract UMI from read           full repeat reads don't need to have LBC and UMI
        if 'null' not in LBC:       
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

            if not LBC_end and not double:     #full repeat reads don't need LBC (it is cut off)
                NonLBCReads+=1
                if contaminations is not None:                                                                                                                            ###### contaminations
                    for contamination in contaminationlist:
                        if contamination in L:
                            contamination_dict[contamination] += 1
                continue
                

            if UMI_status =='on':
                if LBC_end:
                    UMI_present=1
                    umi=L[LBC_end:LBC_end+UMI_length]
                else:
                    UMI_present=0
                    if not double:          #full repeat reads don't need UMI (it is cut off)
                        if contaminations is not None:                                                                                                                            ###### contaminations
                            for contamination in contaminationlist:
                                if contamination in L:
                                    contamination_dict[contamination] += 1
                        continue
            
                if UMI_present == 1 and len(umi) != UMI_length:  
                    NonUMIReads+=1
                    UMI_present=0

            else:
                UMI_present=0                    
            
        if double:
            SinglefullRepeatReads += 1
            I.write(readName+'\n'+L+'\n')
            spacerProximal = ''
            spacerProximal = editSpacer_WFA(tempSpacerProximal,firstRepeat,secondRepeat,minSpacer,maxSpacer, dr1_aligner_dict, dr2_aligner_dict)
            spacerDistal = ''
            tempSpacerDistal = L[(full_end-len(firstRepeat)-8):] 
            if len(tempSpacerDistal) >= minReadLength:
                    spacerDistal = editSpacer_WFA(tempSpacerDistal,firstRepeat,secondRepeat,minSpacer,maxSpacer, dr1_aligner_dict, dr2_aligner_dict) 
                
            if spacerDistal and spacerProximal:
                spacerReads+=1
                spacerReadsDoubleBoth+=1
                doubleSpacer= spacerProximal+"_"+spacerDistal
                    
                # store spacers in dict, this will force uniqueness                       
                if doubleSpacer not in allDistalSpacers:
                    if doubleSpacer not in D:
                        D[doubleSpacer]=[readName+'_doubleAcquisitions_both',0]
                    D[doubleSpacer][1]+=1
                    allDistalSpacers.append(spacerDistal)
                    # J.write(readName+'_doubleAcquisitions_both_distal'+'\n'+spacerDistal+'\n')
                    # J.write(readName+'_doubleAcquisitions_both_proximal'+'\n'+spacerProximal+'\n')
                K.write(readName+'_doubleAcquisitions_both_distal'+'\n'+spacerDistal+'\n')
                K.write(readName+'_doubleAcquisitions_both_proximal'+'\n'+spacerProximal+'\n')
                if UMI_present == 1: 
                    # if umi and spacer are both new #
                    if umi not in UMISpacers and doubleSpacer not in UniqueSpacersUMI: # Check if UMI is currently in the dictionary, if not, add it
                        UniqueUMIs+=1 
                        # UMIcounts[umi]=1 # new UMI
                        UMISpacers[umi]=[[doubleSpacer, 1]] # add umi with spacer to UMI dictionary
                        UniqueSpacersUMI[doubleSpacer]= [umi] # add spacer with umi to Spacer dictionary
                        UniqueSpacersUMIReadNames[doubleSpacer]=[readName+'_doubleAcquisitions_both']
                    # if spacer is in the spacer dictionary, but umi is new #
                    elif umi not in UMISpacers:
                        UniqueUMIs+=1 
                        # UMIcounts[umi]=1
                        UMISpacers[umi]=[[doubleSpacer, 1]]
                        UniqueSpacersUMI[doubleSpacer].append(umi)
                        UniqueSpacersUMIReadNames[doubleSpacer].append(readName+'_doubleAcquisitions_both')
                    # if UMI is currently in the umi dictionary, check if corresponding spacers have very similar sequences
                    elif doubleSpacer not in UniqueSpacersUMI:
                        sameSpacerSignal=0 # signal to turn on if spacer matches
                        # check the umi dictionary for the spacer or closely similar (sub)sequence
                        for i in UMISpacers[umi]:
                            if '_' not in i[0] and len(i[0])>=len(doubleSpacer):   
                                continue
                            if '_' in i[0] and len(i[0])>=len(doubleSpacer): 
                                aligner = WavefrontAligner(doubleSpacer, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                                match = aligner(i[0], clip_cigar=False)
                                if match.score > -30:
                                    sameSpacerSignal=1
                                    i[1]+=1
                                    break
                            else:
                                aligner = WavefrontAligner(i[0], gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                                match = aligner(doubleSpacer, clip_cigar=False)
                                if match.score > -30:                                                   # MAYBE CHECK THESE SCORES
                                    sameSpacerSignal=1
                                    del UniqueSpacersUMIReadNames[i[0]][UniqueSpacersUMI[i[0]].index(umi)]    
                                    if not UniqueSpacersUMIReadNames[i[0]]:
                                        del UniqueSpacersUMIReadNames[i[0]]
                                    del UniqueSpacersUMI[i[0]][UniqueSpacersUMI[i[0]].index(umi)]
                                    if not UniqueSpacersUMI[i[0]]:
                                        del UniqueSpacersUMI[i[0]]
                                    UniqueSpacersUMI[doubleSpacer]=[umi]
                                    UniqueSpacersUMIReadNames[doubleSpacer]=[readName+'_doubleAcquisitions_both']
                                    i[0]=doubleSpacer
                                    i[1]+=1
                                    break
                        if sameSpacerSignal == 0:
                        # if the spacer is new
                            # UMIcounts[umi]+=1
                            UMISpacers[umi].append([doubleSpacer, 1])
                            UniqueSpacersUMI[doubleSpacer]=[umi]
                            UniqueSpacersUMIReadNames[doubleSpacer]=[readName+'_doubleAcquisitions_both']
            # process reads with only distal spacer
            elif spacerDistal:
                spacerReads+=1
                spacerReadsDoubleOne+=1
                if spacerDistal not in allDistalSpacers:
                    if spacerDistal not in D:
                        D[spacerDistal]=[readName+'_doubleAcquisitions_distal',0]
                    D[spacerDistal][1]+=1
                K.write(readName+'_doubleAcquisitions_distal'+'\n'+spacerDistal+'\n')
                if UMI_present == 1:
                    if umi not in UMISpacers and spacerDistal not in UniqueSpacersUMI: # Check if UMI is currently in the dictionary, if not, add it
                        UniqueUMIs+=1
                        # UMIcounts[umi]=1
                        UMISpacers[umi]=[[spacerDistal, 1]]
                        UniqueSpacersUMI[spacerDistal]= [umi]
                        UniqueSpacersUMIReadNames[spacerDistal]=[readName+'_doubleAcquisitions_distal']
                    elif umi not in UMISpacers:
                        UniqueUMIs+=1
                        # UMIcounts[umi]=1
                        UMISpacers[umi]=[[spacerDistal, 1]]
                        UniqueSpacersUMI[spacerDistal].append(umi)
                        UniqueSpacersUMIReadNames[spacerDistal].append(readName+'_doubleAcquisitions_distal')
                    elif spacerDistal not in UniqueSpacersUMI:
                    # if UMI is currently in the dictionary, check if corresponding spacers are the same (with some leeway for errors)
                        sameSpacerSignal=0
                        for i in UMISpacers[umi]:
                            if len(i[0])>=len(spacerDistal):
                                aligner = WavefrontAligner(spacerDistal, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                                match = aligner(i[0], clip_cigar=False)
                                if match.score > -30:
                                    if 'double' not in UniqueSpacersUMIReadNames[i[0]][UniqueSpacersUMI[i[0]].index(umi)]:
                                        UniqueSpacersUMIReadNames[i[0]][UniqueSpacersUMI[i[0]].index(umi)] = [readName+'_doubleAcquisitions_distal']
                                    sameSpacerSignal=1
                                    i[1]+=1
                                    break
                            else:
                                aligner = WavefrontAligner(i[0], gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                                match = aligner(spacerDistal, clip_cigar=False)
                                if match.score > -30:
                                    sameSpacerSignal=1
                                    del UniqueSpacersUMIReadNames[i[0]][UniqueSpacersUMI[i[0]].index(umi)] 
                                    if not UniqueSpacersUMIReadNames[i[0]]:
                                        del UniqueSpacersUMIReadNames[i[0]]  
                                    del UniqueSpacersUMI[i[0]][UniqueSpacersUMI[i[0]].index(umi)]
                                    if not UniqueSpacersUMI[i[0]]:
                                        del UniqueSpacersUMI[i[0]]
                                    UniqueSpacersUMI[spacerDistal]=[umi]
                                    UniqueSpacersUMIReadNames[spacerDistal]=[readName+'_doubleAcquisitions_distal']
                                    i[0]=spacerDistal
                                    i[1]+=1
                                    break
                        if sameSpacerSignal == 0:
                            # UMIcounts[umi]+=1
                            UMISpacers[umi].append([spacerDistal, 1])
                            UniqueSpacersUMI[spacerDistal]=[umi]
                            UniqueSpacersUMIReadNames[spacerDistal]=[readName+'_doubleAcquisitions_distal']
            # process reads with only proximal spacer
            elif spacerProximal:
                spacerReads+=1
                spacerReadsDoubleOne+=1
                if spacerProximal not in D:
                    D[spacerProximal]=[readName+'_doubleAcquisitions_proximal',0]
                D[spacerProximal][1]+=1
                K.write(readName+'_doubleAcquisitions_proximal'+'\n'+spacerProximal+'\n')
                if UMI_present == 1:
                    if umi not in UMISpacers and spacerProximal not in UniqueSpacersUMI: # Check if UMI is currently in the dictionary, if not, add it
                        UniqueUMIs+=1
                        # UMIcounts[umi]=1
                        UMISpacers[umi]=[[spacerProximal, 1]]
                        UniqueSpacersUMI[spacerProximal]= [umi]
                        UniqueSpacersUMIReadNames[spacerProximal]=[readName+'_doubleAcquisitions_proximal']
                    elif umi not in UMISpacers:
                        UniqueUMIs+=1
                        # UMIcounts[umi]=1
                        UMISpacers[umi]=[[spacerProximal, 1]]
                        UniqueSpacersUMI[spacerProximal].append(umi)
                        UniqueSpacersUMIReadNames[spacerProximal].append(readName+'_doubleAcquisitions_proximal')
                    elif spacerProximal not in UniqueSpacersUMI:
                    # if UMI is currently in the dictionary, check if corresponding spacers are the same (with some leeway for errors)
                        sameSpacerSignal=0
                        for i in UMISpacers[umi]:
                            if "_" in i[0]:
                                continue # if spacer proximal is part of earlier double/multiple acquisition, it would still count as unique 
                            # since no distal spacer found and proximal spacer can't be amplified singly without distal
                            if len(i[0])>len(spacerProximal): 
                                aligner = WavefrontAligner(spacerProximal, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                                match = aligner(i[0], clip_cigar=False)
                                if match.score > -30:
                                    if 'double' not in UniqueSpacersUMIReadNames[i[0]][UniqueSpacersUMI[i[0]].index(umi)]:
                                        UniqueSpacersUMIReadNames[i[0]][UniqueSpacersUMI[i[0]].index(umi)] = [readName+'_doubleAcquisitions_distal']
                                    sameSpacerSignal=1
                                    i[1]+=1
                                    break
                            else:
                                aligner = WavefrontAligner(i[0], gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                                match = aligner(spacerProximal, clip_cigar=False)
                                if match.score > -30:
                                    sameSpacerSignal=1
                                    del UniqueSpacersUMIReadNames[i[0]][UniqueSpacersUMI[i[0]].index(umi)]
                                    if not UniqueSpacersUMIReadNames[i[0]]:
                                        del UniqueSpacersUMIReadNames[i[0]]    
                                    del UniqueSpacersUMI[i[0]][UniqueSpacersUMI[i[0]].index(umi)]
                                    if not UniqueSpacersUMI[i[0]]:
                                        del UniqueSpacersUMI[i[0]]
                                    UniqueSpacersUMI[spacerProximal]=[umi]
                                    UniqueSpacersUMIReadNames[spacerProximal]=[readName+'_doubleAcquisitions_proximal']
                                    i[0]=spacerProximal
                                    i[1]+=1
                                    break
                        if sameSpacerSignal == 0:
                            # UMIcounts[umi]+=1
                            UniqueSpacersUMI[spacerProximal]=[umi]
                            UMISpacers[umi].append([spacerProximal, 1])
                            UniqueSpacersUMIReadNames[spacerProximal]=[readName+'_doubleAcquisitions_proximal']
            else:
                if contaminations is not None:                                                                                                                            ###### contaminations
                    for contamination in contaminationlist:
                        if contamination in L:
                            contamination_dict[contamination] += 1
                continue

       # no double acquisition
        else:
            spacer = editSpacer_WFA(L,firstRepeat,secondRepeat,minSpacer,maxSpacer, dr1_aligner_dict, dr2_aligner_dict)
            if spacer == '':
                if contaminations is not None:                                                                                                                         ###### contaminations
                    for contamination in contaminationlist:
                        if contamination in L:
                            contamination_dict[contamination] += 1
                continue
                                                                                                 
            spacerReads+=1

            if spacer not in D:
                D[spacer]=[readName,0]
            D[spacer][1]+=1

            if UMI_present == 1:
                if umi not in UMISpacers and spacer not in UniqueSpacersUMI:
                    UniqueUMIs+=1
                    # UMIcounts[umi]=1
                    UMISpacers[umi]=[[spacer, 1]]
                    UniqueSpacersUMI[spacer]= [umi]
                    UniqueSpacersUMIReadNames[spacer]=[readName]
                
                elif umi not in UMISpacers:
                    UniqueUMIs+=1
                    # UMIcounts[umi]=1
                    UMISpacers[umi]=[[spacer, 1]]
                    UniqueSpacersUMI[spacer].append(umi)
                    UniqueSpacersUMIReadNames[spacer].append(readName)
                
                elif spacer not in UniqueSpacersUMI:
                    sameSpacerSignal=0
                    for i in UMISpacers[umi]:
                        if len(spacer)<=len(i[0]):
                            aligner = WavefrontAligner(spacer, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                            match = aligner(i[0], clip_cigar=False)
                            if match.score > -30:
                                sameSpacerSignal=1
                                i[1]+=1
                                break
                        else:
                            aligner = WavefrontAligner(i[0], gap_opening = 4, gap_opening2 = 24, gap_extension = 2, span="end-to-end")
                            match = aligner(spacer, clip_cigar=False)
                            if match.score > -30:
                                sameSpacerSignal=1
                                del UniqueSpacersUMIReadNames[i[0]][UniqueSpacersUMI[i[0]].index(umi)]
                                if not UniqueSpacersUMIReadNames[i[0]]:
                                    del UniqueSpacersUMIReadNames[i[0]]
                                del UniqueSpacersUMI[i[0]][UniqueSpacersUMI[i[0]].index(umi)]
                                if not UniqueSpacersUMI[i[0]]:
                                    del UniqueSpacersUMI[i[0]]
                                UniqueSpacersUMI[spacer]=[umi]
                                UniqueSpacersUMIReadNames[spacer]=[readName]
                                i[0]=spacer
                                i[1]+=1
                                break

                    if sameSpacerSignal == 0:
                        # UMIcounts[umi]+=1
                        UMISpacers[umi].append([spacer, 1])
                        UniqueSpacersUMI[spacer]=[umi]
                        UniqueSpacersUMIReadNames[spacer]=[readName]

            K.write(readName+'\n'+spacer+'\n')

    # iterate through spacers and print to file
    seqlist=sorted(D.keys())    
    for S in seqlist:
        if '_doubleAcquisitions_both' in D[S][0]:
            UniqueDoubleAcquisitions+=1
        # elif 'multiple' in D[S][0]:
        #     UniqueMultipleAcquisitions+=1
        else:
            UniqueSingleAcquisitions+=1
        if '_' not in S and len(S)>0:
            G.write(D[S][0]+'_rep_'+str(D[S][1])+'\n'+S+'\n')
            count = Counter(S)
            gc = (count['G']+count['C'])*100/len(S)
            if infoFile:
                GC.write(str(S)+'\t'+str(len(S))+'\t'+str(round(gc,2))+'\t'+str(D[S][1])+'\n')

        else:
            for index, spcr in enumerate(S.split('_')):
                if len(spcr)>0:
                    G.write(D[S][0]+'_rep_'+str(D[S][1])+'_spacerPosition_'+str(index)+'\n'+spcr+'\n')
                    count = Counter(spcr)
                    gc = (count['G']+count['C'])*100/len(spcr)
                    if infoFile:
                        GC.write(str(spcr)+'\t'+str(len(spcr))+'\t'+str(gc)+'\n')


    if UMI_status =='on':

        for sq in UniqueSpacersUMIReadNames:
            for readname in UniqueSpacersUMIReadNames[sq]:
                UniqueAcquisitionsUMI+=1
                if '_doubleAcquisitions_both' in readname:
                    UniqueDoubleAcquisitionsUMI+=1
                # elif 'multiple' in readname:
                #     UniqueMultipleAcquisitionsUMI+=1
                else:
                    UniqueSingleAcquisitionsUMI+=1   
                if '_' not in sq and len(sq)>0:
                    U.write(str(readname)+'\n'+str(sq)+'\n')
                    # count = Counter(sq)
                    # gc = (count['G']+count['C'])*100/len(sq)
                    # if infoFile:
                    #     UGC.write(str(sq)+'\t'+str(len(sq))+'\t'+str(gc)+'\n')
                else:
                    for index, spcr in enumerate(sq.split('_')):
                        if len(spcr)>0:
                            U.write(readname+'_spacerPosition_'+str(index)+'\n'+spcr+'\n')
                            # count = Counter(spcr)
                            # gc = (count['G']+count['C'])*100/len(spcr)
                            # if infoFile:
                            #     UGC.write(str(spcr)+'\t'+str(len(spcr))+'\t'+str(gc)+'\n')

        for umi in UMISpacers:
            for spcr in UMISpacers[umi]:
                count = Counter(spcr[0])
                gc = (count['G']+count['C'])*100/len(spcr[0])
                if infoFile:
                    if '_' not in spcr[0]:
                        UGC.write(str(spcr[0])+'\t'+str(umi)+'\t'+str(len(spcr[0]))+'\t'+str(round(gc,2))+'\t'+str(spcr[1]) +'\n')
                    else:
                        UGC.write(str(spcr[0])+'\t'+str(umi)+'\t'+str(len(spcr[0])-1)+'\t'+str(round(gc,2))+'\t'+str(spcr[1]) +'\n')

    if not D:
        os.system(str("echo 'No spacers to map'"))
    os.system(str('echo '+"'"+'*'+'\t'+' sampleName'+'\t'+' rawReads'+'\t'+'identifiedSpacers'+'\t'+'uniqueSpacers'+'\t'+'spacersWithFullDR'+'\t'+'doubleAcquisitions'+'\t'+'uniqueSingleAcquisitions'+'\t'+'uniqueDoubleAcquisitions'+"'"))
    # os.system(str('echo '+"'"+'@'+'\t'+str(outPath+outName)+'\t'+str(rawReads)+'\t'+str(spacerReads)+'\t'+str(len(D.keys()))+'\t'+str(SinglefullRepeatReads)+'\t'+str(spacerReadsDoubleBoth)+'\t'+str(MultiplefullRepeatReads)+'\t'+str(spacerReadsMultiComplete)+'\t'+str(UniqueSingleAcquisitions)+'\t'+str(UniqueDoubleAcquisitions)+'\t'+str(UniqueMultipleAcquisitions)+"'"))
    # SS.write(str(outName)+'\t'+str(rawReads)+'\t'+str(spacerReads)+'\t'+str(len(D.keys()))+'\t'+str(SinglefullRepeatReads)+'\t'+str(spacerReadsDoubleBoth)+'\t'+str(MultiplefullRepeatReads)+'\t'+str(spacerReadsMultiComplete)+'\t'+str(UniqueSingleAcquisitions)+'\t'+str(UniqueDoubleAcquisitions)+'\t'+str(UniqueMultipleAcquisitions)+'\n')
    os.system(str('echo '+"'"+'@'+'\t'+str(outPath+outName)+'\t'+str(rawReads)+'\t'+str(spacerReads)+'\t'+str(len(D.keys()))+'\t'+str(SinglefullRepeatReads)+'\t'+str(spacerReadsDoubleBoth)+'\t'+str(UniqueSingleAcquisitions)+'\t'+str(UniqueDoubleAcquisitions)+"'"))
    SS.write(str(outName)+'\t'+str(rawReads)+'\t'+str(NonLBCReads)+'\t'+str(spacerReads)+'\t'+str(len(D.keys()))+'\t'+str(SinglefullRepeatReads)+'\t'+str(UniqueSingleAcquisitions)+'\t'+str(UniqueDoubleAcquisitions)+'\n')
    
    if contaminations is not None:
        contamination_stats = [str(outName)]                                                                                                                                   ###### contaminations
        for contamination, count in contamination_dict.items():
            if rawReads > 0:
                percentage = (count / rawReads) * 100
            else:
                percentage = 0 

            contamination_stats.append(f'{count}\t{percentage:.2f}')  # Format percentage to 2 decimal places
        row_string = '\t'.join(contamination_stats) + '\n'
        CS.write(row_string)

    if UMI_status=='on':
        # US.write(str(outName)+'\t'+str(rawReads)+'\t'+str(UniqueUMIs)+'\t'+str(NonUMIReads)+'\t'+str(ConstantRegionNonUMIReads)+'\t'+str(UniqueAcquisitionsUMI)+'\t'+str(UniqueSingleAcquisitionsUMI)+'\t'+str(UniqueDoubleAcquisitionsUMI)+'\t'+str(UniqueMultipleAcquisitionsUMI)+'\n')
        US.write(str(outName)+'\t'+str(rawReads)+'\t'+str(NonLBCReads)+'\t'+str(UniqueUMIs)+'\t'+str(NonUMIReads)+'\t'+str(UniqueAcquisitionsUMI)+'\t'+str(UniqueSingleAcquisitionsUMI)+'\t'+str(UniqueDoubleAcquisitionsUMI)+'\n')

    
    # print useful summary stats to stdout

    F.close()
    G.close()
    I.close()
    # J.close()
    K.close()
    # MC.close()
    # M.close()
    SS.close()
    if contaminations is not None:
        CS.close()                                                               ###### contaminations
    if infoFile:
        GC.close()
    if UMI_status =='on':
        U.close()
        # UR.close()
        if infoFile:
            UGC.close()