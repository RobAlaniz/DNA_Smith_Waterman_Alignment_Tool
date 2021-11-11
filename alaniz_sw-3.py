# imports
import pandas as pd
import numpy as np
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.pairwise2 import format_alignment
from Bio import pairwise2 as pw2
import argparse
import sys
import os



# arguments that are command line ready


parser = argparse.ArgumentParser(description = "Problem Set 2: Smith-Waterman Local Alignment.")
parser.add_argument('-g', type=int, required=False)
parser.add_argument('-m', type=str, required=False)
parser.add_argument("-s1", metavar="s1", action='store', type=str, help="Sequence 1")
parser.add_argument("-s2", metavar="s2", action='store', type=str, help="Sequence 2")
parser.add_argument("-o", action="store", type=str, help="Output File")

args = parser.parse_args()

output_file = open(args.o, 'a+')


# I am initializing a bunch of variables here

gap = 0

seq1_string = ""
seq2_string = ""

matrix_string = ""
matrix_string = str(args.m)

m_string = str(matrix_string)

blosum_matrix_file = ""


print("\n")

# This section of code is to verify the users input on command line.
# If the user does not use the correct notation then the program will
# run the exit() call and end the program

if m_string == "blosum62" or m_string == "BLOSUM62" or m_string == "Blosum62" or m_string == "blosum62.txt" or m_string == "Blosum62.txt" or m_string == "BLOSUM62.txt":
    blosum_matrix_file = "BLOSUM62.txt"
elif m_string == "blosum50" or m_string == "BLOSUM50" or m_string == "Blosum50" or m_string == "blosum50.txt" or m_string == "Blosum50.txt" or m_string == "BLOSUM50.txt":
    blosum_matrix_file = "BLOSUM50.txt"
else:
    
    user_input = input("Choose a Matrix:\nBLOSUM50: Press 1 then press enter.\n\nBLOSUM62: Press 2 then press enter.\n\nChoice: ")
    if user_input == str(1):
        blosum_matrix_file = "BLOSUM50.txt"
        print("\nComputing using BLOSUM50 matrix.")
    elif user_input == str(2):
        blosum_matrix_file = "BLOSUM62.txt"
        print("\nComputing using BLOSUM62 matrix.")
    else:
        print("Invalid input, the program will end.")
        exit()



# Args cannot be directly casted to another data type and need to go through a transitional faze
# in order to create a data type that is usable within the program

if args.g == None:
    gap = 8
else:
    gap_string = str(args.g)
    gap = int(gap_string)







# The fasta files are read in here and stored as seq_1 and seq_2 respectively.

def read_file(filename):
    file = open(filename, "r")
    seq_array = file.readlines()
    seq_output = []
    file.close()
    for i in seq_array[1:]:
        seq_output += list(i.strip())
    return seq_output

#This is where I created string files that correspond to both of the sequences
print("\n")
if args.s1:
    seq_1 = read_file(args.s1)
if args.s2:
    seq_2 = read_file(args.s2)
for i in seq_1:
    seq1_string = seq1_string + i


for i in seq_2:
    seq2_string = seq2_string + i











#Opening the BLOSUM File and iterating through the code Smith-Waterson Algorithm.
#Either BLOSUM50 or BLOSUM62 depending on what was either entered or chosen from the menu

f = open(blosum_matrix_file, "r")
blosum_file_data = f.readlines()[6:]

array = []
for a in blosum_file_data:
    b = a.strip().split()
    array.append(b)

computed_array = []
for c in array[1:]:
    computed_array.append(c[1:])

blosum_dataframe = pd.DataFrame(computed_array, columns = array[0], index = array[0]).apply(pd.to_numeric)

# This sets all the scores in the matrix to 0
new_matrix = np.zeros((len(seq_1)+1,len(seq_2)+1))

#this function iterates through the matrix scoring as it goes through the individual cells
def score(new_matrix, seq_1, seq_2, i, j):
    result = [0, new_matrix[i-1][j]-gap, new_matrix[i][j-1]-gap, new_matrix[i-1][j-1] + blosum_dataframe[seq_1[i-1]][seq_2[j-1]]]

    return max(result)


for i in range(1,len(seq_1)+1):
    for j in range(1,len(seq_2)+1):

        new_matrix[i][j] = score(new_matrix, seq_1, seq_2, i, j)

# While the length and width of the matrix are less than 100
# the Maximum Individual score is calculated and located
if len(seq_1) <= 100 and len(seq_2) <= 100:
    print("\nTHE RESULTING SCORE MATRIX:\n")
    print(new_matrix)
    print("Score Matrix:\n", file=output_file)
    print(new_matrix, file=output_file)

(max_i,max_j) = np.unravel_index(np.argmax(new_matrix), new_matrix.shape)
print("\n")
print("MAXIMUM INDIVIDUAL SCORE: ")
max_score = np.max(new_matrix)
print(max_score, "\n")


# This is where I used BioPython to read in a alignment methods
# I was also able to set the match, mismatch, and gap scores

for a in pw2.align.localxx(seq1_string, seq2_string):
    print(format_alignment(*a, full_sequences=True))
    
aligner = Align.PairwiseAligner()
print("\n")
print("\nUsing " + blosum_matrix_file +  " matrix file.\n")


# This is where I calculate the Alignments for the sequence comparison
# depending on the User's Matrix Choice, the alignments will be scored
# and the Maximum Individual Score within the chosen matrices is found.
# I calculate the Similarity Score, Alignment, and percent similarity in
# this section of the code

if blosum_matrix_file == "BLOSUM62.txt":
    
    aligner = Align.PairwiseAligner()
    aligner.gap_score = (-1 * gap)
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    alignments = aligner.align(seq1_string, seq2_string)
    alignments = list(alignments)
    alignment = alignments[0]
    print("Match Score = %.f" % alignment.score)
    local_align = pw2.align.localxx(seq1_string, seq2_string)
    seq_length = min(len(seq1_string), len(seq2_string))
    matches = local_align[0][2]
    percent_match = (matches / seq_length) * 100
    print(percent_match, "%")

elif blosum_matrix_file == "BLOSUM50.txt":
    aligner = Align.PairwiseAligner()
    aligner.gap_score = gap
    aligner.match_score = 2
    aligner.mismatch_score = -1
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM50")
    alignments = aligner.align(seq1_string, seq2_string)
    alignments = list(alignments)
    alignment = alignments[0]
    print("Match Score = %.f" % alignment.score)
    local_align = pw2.align.localxx(seq1_string, seq2_string)
    seq_length = min(len(seq1_string), len(seq2_string))
    matches = local_align[0][2]
    percent_match = (matches / seq_length) * 100
    print(percent_match, "%")
    
print("\n")


output_file.close()
