import csv
import sys, os
import numpy as np
import math, statistics

# Define the directory containing the fasta files and the sequence to find
fasta_directory = "/cluster/scratch/hugifl/quick_test/"
files = os.listdir(fasta_directory)
dr2 = 'GTCGTACTTT'

# Initialize lists to store start and end positions of matches
start_positions = []
end_positions = []

# Loop through each file in the directory
for file in files:
    if file.endswith(".fasta"):  # Ensure we are reading fasta files
        with open(fasta_directory + file, 'r') as f:
            # Read the file line by line
            for line in f:
                if not line.startswith('>'):  # Ignore header lines
                    # Find all occurrences of dr2 in the current sequence
                    start = 0
                    while start != -1:
                        start = line.find(dr2, start)
                        if start != -1:
                            end = start + len(dr2)
                            start_positions.append(start)
                            end_positions.append(end)
                            start += 1  # Move start to the next position for searching

# Calculate mean and standard deviation for start and end positions
mean_start = np.mean(start_positions)
std_dev_start = np.std(start_positions)
mean_end = np.mean(end_positions)
std_dev_end = np.std(end_positions)

# Print the results
print(f"Mean start position: {mean_start}, Standard Deviation: {std_dev_start}")
print(f"Mean end position: {mean_end}, Standard Deviation: {std_dev_end}")

