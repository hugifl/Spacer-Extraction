# Oct, 2023  Florian Hugi 
# Function to construct WFA2 aligner objects based on a minimal and maximal length of the (part of the) read where a spacer is searched. 

from pywfa import WavefrontAligner

def construct_aligners(pattern, min_text_length, max_text_length):
    aligner_dict = {}
    for length in range(min_text_length, max_text_length + 1, 10):
        aligner = WavefrontAligner(pattern, gap_opening = 4, gap_opening2 = 24, gap_extension = 2, text_end_free = length, text_begin_free = length)
        aligner_dict[length] = aligner

    return aligner_dict