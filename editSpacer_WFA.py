# Oct, 2023  Florian Hugi 
# Function to extract a spacer from a predefined part of a read using WFA2 aligner objects. 
# Adapted from Tanmay Tanna

import collections

# define a custom data type for found matches
Match = collections.namedtuple('Match', ['start', 'end', 'dist'])

def editSpacer_WFA(read,firstRepeat,secondRepeat,minSpacer,maxSpacer, dr1_aligner_dict, dr2_aligner_dict): 
    findDR1 = read.find(firstRepeat)
    firstMatch=''
    spacerStart=''
    if findDR1 ==-1:
        dr1_aligner = dr1_aligner_dict[(len(read) // 10) * 10]
        firstMatch = dr1_aligner(read, clip_cigar=False)
        if (firstMatch.text_end - firstMatch.text_start)>=6 and firstMatch.score > -8:
            spacerStart = firstMatch.text_end 
        
    elif findDR1 !=-1:
        
        firstMatch = [Match(findDR1, findDR1 + len(firstRepeat),0)]
        spacerStart = firstMatch[0].end 
    
    if not spacerStart:  
        return ''   # No alignment. Return empty string.
        
    findDR2 = read.find(secondRepeat)
    secondMatch=''
    spacerEnd = ''
    if findDR2 ==-1:
        dr2_aligner = dr2_aligner_dict[(len(read) // 10) * 10]
        secondMatch = dr2_aligner(read, clip_cigar=False)
        if (secondMatch.text_end - secondMatch.text_start)>=6 and secondMatch.score > -8:
            spacerEnd = secondMatch.text_start 

    elif findDR2 != -1: 
        secondMatch = [Match(findDR2, findDR2 + len(secondRepeat),0)]
        spacerEnd = secondMatch[0].start 
    if not spacerEnd: 
        return ''

    spacer = read[spacerStart:spacerEnd]

    # no spacer if out of bounds
    if len(spacer) > maxSpacer: return ''
    if len(spacer) < minSpacer: return ''
    
    return spacer