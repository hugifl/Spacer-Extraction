# This function searches full repeats in the read and checks if a full repeat (double acquisition) is present.
# Because in some cases, full repeats can appear early in the read. In this case, the read is truncated and the function again searches for a full repeat.


def check_double(L, fullRepeat, drfull_aligner_dict, minReadLength, secondRepeat):
    double = False
    skip_read = False
    tempSpacerProximal = ''
    full_end = ''
    findfull = L.find(fullRepeat)
    if findfull ==-1:
        aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
        drfull_match = aligner_drfull(L, clip_cigar=False)
        full_start = drfull_match.text_start
        full_end = drfull_match.text_end
        double = drfull_match.score > -12 and (drfull_match.text_end - drfull_match.text_start) >= 25
    elif findfull !=-1:
        full_start = findfull
        full_end = findfull + len(fullRepeat)
        double = True
    
    if double:  # if full repeat is present
        # split read into (leader) proximal and (leader) distal and independently run the editSpacer code to extract proximal/distal spacers
        tempSpacerProximal = L[:(full_start+len(secondRepeat))+8] #8 more nucleotides are added to be sure to have the full secondRepeat in the string
        if len(tempSpacerProximal) < minReadLength:  
            L = L[len(tempSpacerProximal) - 5:]
            trim_L = True
            double = False
            if len(L) >= minReadLength:
                aligner_drfull = drfull_aligner_dict[(len(L) // 10) * 10]
                drfull_match = aligner_drfull(L, clip_cigar=False)
                full_start = drfull_match.text_start
                full_end = drfull_match.text_end
                double = drfull_match.score > -14 and (drfull_match.text_end - drfull_match.text_start) >= 25
                tempSpacerProximal = L[:(full_start+len(secondRepeat))+8]
                if len(tempSpacerProximal) < minReadLength:
                    double = False
            else:
                skip_read = True
    
    return skip_read, double, tempSpacerProximal, full_end, L

