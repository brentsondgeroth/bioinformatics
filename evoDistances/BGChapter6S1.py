'''
Name: Brent Gaither Kevin Portland Lei Guo
Description: This program uses the Jukes-Cantor model to estimate evolutionary distances.
Due Date: 2/10/15

'''
import math

'''
readIn- Reads in the files
'''

def readIn(infile):
    sequence = ""
    infile.readline()  # bypass > header line
    for line in infile:
        line = line.replace('\n', '')
        sequence = sequence + line
    return sequence.upper()

'''
compareSequences- reads in the sequences and compares to find matches and length of sequences.
'''

def compareSequences(sequence1,sequence2):
	lenCtr = 0
	difCtr = 0
	for char1, char2 in zip(sequence1, sequence2):
		if char1 != "-" and char2 != "-": #Does not count unmatched nucleotides. 
			lenCtr += 1
			if char1!= char2:
				difCtr += 1
	return lenCtr,difCtr
'''
JukesCantor- Uses the differences and lengths of nucleotide sequences to 
find evolutionary distances.
'''
def jukesCantor(difCtr,lenCtr):
	difCtr = difCtr / lenCtr #converts to fraction of substitutions 
	jukes = ((-3/4) * math.log(1-(4/3*difCtr))) 
	return jukes

'''
Main
'''

infile1 = open('hippo.txt', 'r')
infile2 = open('giraffe.txt', 'r')
outfile = open('ch6sk1out.txt', 'w')

sequence1 = readIn(infile1)
sequence2 = readIn(infile2)
infile1.close()
infile2.close()

lenCtr, difCtr = compareSequences(sequence1,sequence2)
jukes = jukesCantor(difCtr,lenCtr)

print("This program calculates the evolutionary distance between two sequences using the Jukes-Cantor model.")
print(jukes)
outfile.write("Jukes-Cantor model distance \n")
outfile.write(str(jukes))
outfile.close