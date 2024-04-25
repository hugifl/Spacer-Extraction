# Nov 15, 2022
# Tanmay Tanna and Wiona Glaenzer

import fuzzysearch
import collections

# define a custom data type for found matches
Match = collections.namedtuple('Match', ['start', 'end', 'dist'])

def editSpacer(read,firstRepeat,firstExpect,secondRepeat,secondExpect,DR1Mismatch,DR2Mismatch,minSpacer,maxSpacer,firstRangeToLook,secondRangeToLook):
    s=''
    # Find the first repeat using string search, if not found, use fuzzy string matching

    s = read[firstExpect:(firstExpect+firstRangeToLook+len(firstRepeat))]
    #print("first" +str(len(s))) #--------------------------------------- print -------------------------
    findDR1 = s.find(firstRepeat)
    firstMatch=''
    if findDR1 ==-1 and DR1Mismatch>0 :
        firstMatch = fuzzysearch.find_near_matches(firstRepeat, s, max_l_dist = DR1Mismatch)
        firstMatch = sorted(firstMatch, key=lambda x: x.dist)
    elif findDR1 !=-1:
        #firstMatch = [[findDR1, findDR1 + len(firstRepeat)]]
        firstMatch = [Match(findDR1, findDR1 + len(firstRepeat),0)]


    if not firstMatch:   return ''   # Too many mismatches. Return empty string.

    s=''
    # Find the second repeat using fuzzy string matching 
   
    s = read[secondExpect:(secondExpect+secondRangeToLook+len(secondRepeat))]
    #print("second" +str(len(s))) #--------------------------------------- print -------------------------
    if len(s)>len(secondRepeat):
        findDR2 = s.find(secondRepeat)
        secondMatch=''
        if findDR2 ==-1 and DR2Mismatch>0:
            secondMatch = fuzzysearch.find_near_matches(secondRepeat, s, max_l_dist = DR2Mismatch)
            secondMatch = sorted(secondMatch, key=lambda x: x.dist)
        elif findDR2 != -1: 
            secondMatch = [Match(findDR2, findDR2 + len(secondRepeat),0)]
        if not secondMatch: return ''
    else:
        return ''
    
       # Too many mismatches or read too short. Return empty string.

    # If both repeats seem to have been found, return spacer
    spacerStart = firstMatch[0].end
    spacerEnd = secondMatch[0].start + secondExpect


    

    # if spacer is too short, looking for other matches for second DR in the read
    if len(secondMatch) > 1 & spacerEnd - spacerStart < minSpacer:
        if findDR2 ==-1:
            i=1
            while i < len(secondMatch) and spacerEnd - spacerStart < minSpacer:
                spacerEnd = secondMatch[i].start + secondExpect
                i+=1
        else:
            secondMatch = fuzzysearch.find_near_matches(secondRepeat, s, max_l_dist = DR2Mismatch)
            secondMatch = sorted(secondMatch, key=lambda x: x.dist)
            i=0
            while i < len(secondMatch) and spacerEnd - spacerStart < minSpacer:
                spacerEnd = secondMatch[i].start + secondExpect
                i+=1

    spacer = read[spacerStart:spacerEnd]
    # DR2seq = read[spacerEnd:spacerEnd+15]
    # DR2.write(DR2seq+'\n')

    # no spacer if out of bounds
    if len(spacer) > maxSpacer: return ''
    if len(spacer) < minSpacer: return ''
    
    return spacer