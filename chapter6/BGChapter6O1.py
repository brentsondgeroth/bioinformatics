'''
Name: Brent Gaither
Description: This program uses the Jukes-Cantor, Tamura and Kimura models 
to estimate evolutionary distances.
Due Date: 2/12/15
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
        line = line.replace('\r', '')
        sequence = sequence + line
    sequence = sequence.upper()
    return sequence

'''
compareSequences- reads in the sequences and compares to find matches, length of sequences,
mutation type and the GC content. 
'''

def compareSequences(alSeq1,alSeq2):
	#Creates and initializes variables 
	lenCtr = float(0) 	
	difCtr = float(0) 	
	gapCtr = float(0) 	
	gcSeq1 = float(0)
	gcSeq2 = float(0)
	transversion = float(0) 	
	transition = float(0)
	gcContent1 = float(0)
	gcContent2 = float(0)
	S = float(0)
	V = float(0)

	mutationType = {"GC":'transversion',"CG":'transversion', "AT":'transversion',"TA":'transversion',
	"AG":'transition',"GA":'transition',"CT":'transition',"TC":'transition',"GT":'transversion',"CA":'transversion',
	"AC":'transversion',"TG":'transversion'} #Creates a dictionary to determine relationship of mutation
	for char1, char2 in zip(alSeq1, alSeq2): #traverses through sequences
		if char1 != "-" and char2 != "-":
			lenCtr += 1
			if char1!= char2:
				difCtr += 1
				key = char1 + char2
				if mutationType[key] == 'transition': #Determines relationship between nucleotide mutation
					transition +=1
				else:
					transversion +=1
			if (char1 == "G" or char1 == "C") and char2 !="-": #Adds up the number of Gs and Cs in sequence
				gcSeq1 += 1
			if (char2 == "G" or char2 == "C") and char2 !="-":
				gcSeq2 += 1

	gcContent1 = (gcSeq1 / lenCtr)
	gcContent2 = (gcSeq2 / lenCtr)
	V = transversion/lenCtr
	S = transition / lenCtr

	return S,V,gcContent1, gcContent2,difCtr,lenCtr

'''
buildMatrix- Builds the matrix for a semi global alignment
'''

def buildMatrix(sequence1,sequence2):

	#sets gap match and mismatch scores
	gapScore = raw_input("Please enter gap score ")
	misMatch = raw_input("Enter mismatch score ")
	match = raw_input("Enter match score ")
	N = len(sequence1)
	M = len(sequence2)
	gapScore = int(gapScore)
	misMatch = int(misMatch)
	match = int(match)
	matrix = [[0 for col in range(M+1)] for row in range(N+1)]#build matrix for semi global alignment
	matrix [0][0] = 0
	for i in range(1,N+1):
		matrix [i][0] = matrix[i-1][0] #No penalty for initial or terminal gaps
		for j in range(1,M+1):
			matrix [0][j] = matrix[0][j-1]
			if (sequence1[i-1] == sequence2[j - 1]):
				score1 = matrix[i - 1][j - 1] + int(match)            
			else:
				score1 = matrix[i - 1][j - 1] + int(misMatch)
			score2 = matrix[i][j - 1] + gapScore
			score3 = matrix[i-1][j] + gapScore
			matrix[i][j] = max(score1, score2, score3)
		dString = buildDirectional(matrix, N, M, gapScore) #Uses buildDirectional to create a directional string 
		alSeq1,alSeq2 = analyzeMatrix(matrix,dString,N,M) #Receives the aligned sequences 
	return alSeq1,alSeq2

'''
buildDirectional- Function that builds a string that provides direction for tracing
a path back through the matrix to build the best match.
'''

def buildDirectional(matrix, N, M, gapScore):
    dString = ''
    currentRow = N
    currentCol = M
    while(currentRow != 0 or currentCol != 0):
        if(currentRow == 0):
            dString = dString + ('H')
            currentCol = currentCol - 1
        elif(currentCol == 0):
            dString = dString + ('V')
            currentRow = currentRow - 1
        elif(matrix[currentRow][currentCol - 1] + int(gapScore) == matrix[currentRow][currentCol]):
            dString = dString + ('H')
            currentCol = currentCol - 1
        elif(matrix[currentRow - 1][currentCol] + int(gapScore) == matrix[currentRow][currentCol]):
            dString = dString + ('V')
            currentRow = currentRow - 1
        else:
            dString = dString + ('D')
            currentRow = currentRow - 1
            currentCol = currentCol - 1
    return dString

'''
analyzeMatrix- Uses the matrix and dString to traverse through and 
place gaps for optimal alignment.
'''

def analyzeMatrix(matrix,dString,N,M):
	dirPos = 0
	alSeq1 = ''
	alSeq2 = ''
	seq1Pos = N -1
	seq2Pos = M -1
	while(dirPos < len(dString)):    
		if(dString[dirPos] == "D"): #Align sequence for match
		    alSeq1 = sequence1[seq1Pos] + alSeq1
		    alSeq2 = sequence2[seq2Pos] + alSeq2
		    seq1Pos = seq1Pos - 1
		    seq2Pos = seq2Pos - 1
		elif(dString[dirPos] == "V"):#Place a space for a gap
		    alSeq1 = sequence1[seq1Pos] + alSeq1
		    alSeq2 = '-' + alSeq2
		    seq1Pos = seq1Pos - 1
		else:
		    alSeq1 = '-' + alSeq1 #Place a space for a gap
		    alSeq2 = sequence2[seq2Pos] + alSeq2
		    seq2Pos = seq2Pos -1
		dirPos = dirPos + 1

	return alSeq1,alSeq2

'''
findJukesCantor- Finds the evolutionary distance between sequences using 
Jukes-Cantor, a one parameter model.
'''

def findJukesCantor(difCtr,lenCtr):
	difCtr = difCtr / lenCtr #converts to fraction of substitutions 
	jukes = ((-3.0/4) * math.log(1-(4.0/3*difCtr)))
	return jukes

'''
findKimura- Finds the evolutionary distance between sequences using 
Kimura, a two parameter model.
'''

def findKimura(S,V,lenCtr):

	kimura = ((1.0/2) * math.log(1/(1-(2*S)-V))) + ((1.0/4) * math.log(1/(1-(2*V))))
	return kimura

'''
findTamura- Finds the evolutionary distance between sequences using 
Tamura, a three parameter model.
'''

def findTamura(S,V,lenCtr,gcContent1,gcContent2):

	C = gcContent1 + gcContent2 - (2 * gcContent1 * gcContent2)

	tamura = -C * math.log(1 - (S/C) - V) - (((1.0/2)* (1-C)) * math.log(1-(2*V)))
	return tamura

'''
Main
'''

infile1 = open('giraffe.txt', 'r')
infile2 = open('whale.txt', 'r')
outfile = open('ch6O1out.txt', 'w')

sequence1 = readIn(infile1)
sequence2 = readIn(infile2)
infile1.close()
infile2.close()
choice = ''

print("This program aligns and calculates the evolutionary distance of two DNA sequences.")

if len(sequence1) != len(sequence2): #aligns sequence if unaligned 
	sequence1,sequence2 = buildMatrix(sequence1,sequence2)

S,V,gcContent1, gcContent2,difCtr,lenCtr = compareSequences(sequence1,sequence2)

while(choice != '4'): 
	choice = raw_input("Please chose a distance calculation. Jukes-Cantor(1) Kimura(2) Tamura(3) exit(4) ")
	if choice == '1':
		jukes = findJukesCantor(difCtr,lenCtr)
		print(jukes)
		outfile.write("Jukes-Cantor evolutionary distance model \n" + str(jukes) + '\n')
	elif choice == '2':
		kimura = findKimura(S,V,lenCtr)
		print(kimura)
		outfile.write("Kimura evolutionary distance model \n" + str(kimura)+ '\n')
	elif choice == '3':
		tamura = findTamura(S,V,lenCtr,gcContent1,gcContent2)
		print(tamura)
		outfile.write("Tamura evolutionary distance model \n" + str(tamura)+ '\n')
outfile.close()
