'''
Name: Brent Gaither
Description: This program uses a training set to create a substitution matrix. Use sequences 
separated by new lines to differentiate sequences. Do not include headers for sequences. 
Place sequences one after the other to be read in as pairs. 
Due Date: 2/5/15
'''
import math

'''
readInTrainingSet- Reads in the files for a training set 
'''

def readInTrainingSet(infile):
	bigSeq1 = []
	bigSeq2 = []
	matched = ''
	with open('ecoliTrainSet.txt') as f: #Opens and closes training set
		for line in f:
			sequence1 = infile.readline()
			sequence2 = infile.readline()
			sequence1 = sequence1.strip().upper()
			sequence2 = sequence2.strip().upper()
			bigSeq1,bigSeq2 = alignSequences(sequence1, sequence2,bigSeq1,bigSeq2)
	return bigSeq1, bigSeq2

'''
alignSequences- reads in the sequences and compares to find matches between sequences.
'''    

def alignSequences(sequence1,sequence2,bigSeq1,bigSeq2):
	seq1List = list(sequence1)
	seq2List = list(sequence2)
	matched = 0
	i = 0
	for char1, char2 in zip(seq1List, seq2List):
	    if char1 == "-" or char2 == "-":
	    	seq1List.pop(i) #delete char that do not have match 
	    	seq2List.pop(i)
	    i += 1
	bigSeq1 = bigSeq1 + seq1List
	bigSeq2 = bigSeq2 + seq2List
	print(bigSeq1)
	print(bigSeq2)
	return bigSeq1, bigSeq2

'''
compareSequences- Counts up and compares the amino acids in the total training set.
'''

def compareSequences(matchMatrix,bigSeq1,bigSeq2,aminoAcidArray):

	values = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
	dictionary = dict(zip(aminoAcidArray, values))
	totalAA = 0
	countDict = {'A':1,'R':1,'N':1,'D':1,'C':1,'Q':1,'E':1,'G':1,'H':1,'I':1,'L':1,'K':1,'M':1,'F':1,'P':1,'S':1,'T':1,'W':1,'Y':1,'V':1}
	for char1, char2 in zip(bigSeq1, bigSeq2):
		matchMatrix[(int(dictionary[char1]))][(int(dictionary[char2]))] += 1
		countDict[char1] += 1 #Adds up number of times a character has been in the sequence
		countDict[char2] += 1
		totalAA += 2
	return countDict,matchMatrix,dictionary, totalAA

'''
createMatrix- Creates a matrix filled with 1s. Then fills the first row and column with 1 letter amino acids.
'''

def createMatrix():
	matrix = [[1 for col in range(21)] for row in range(21)]#Fills matrix with 1's to initialize and uses a pseudo count
	for i in range(0,21):
		matrix[i][0] = aminoAcidArray[i]
		for j in range(0,21):
			matrix[0][j] = aminoAcidArray[j]
	return matrix
'''
createSubMatrix- Creates and computes calculations for the substitution matrix. 
'''
def createSubMatrix(aminoAcidsT,matchMatrix,dictionary,totalAA,countDict):
	qij = 0
	pi = 0
	pj = 0
	eij = 0
	i = 0
	for char1 in aminoAcidsT: #traverse through amino acids and place log value in each
		j = 0
		for char2 in aminoAcidsT:
			qij = ((matchMatrix[(int(dictionary[char1]))][(int(dictionary[char2]))])* 2 / totalAA ) # find qij multiple by 2 to find only matched pairs
			pi =  countDict[char1] /totalAA #number of times char1 is found in sequence 
			pj = countDict[char2] /totalAA
			if char1 == char2: #Used to differentiate if pi and pj are the same
				eij = pi * pj
			else:
				eij = 2 * pi * pj
			subMatrix[(int(dictionary[char1]))][(int(dictionary[char2]))] = math.log(qij/eij,2)
		i+=1
	return subMatrix

'''
main
'''

# Opens files
infile1 = open('ecoliTrainSet.txt', 'r')
outfile = open('ch5O1out.txt', 'w')

aminoAcids  = ("0,A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V") # amino acids used to create matrix
aminoAcidL = ("A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V")#amino acid list used for out file
aminoAcidArray = aminoAcids.split(',')
aminoAcidsT = ('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'); #tuple of amino acids used for dict

matchMatrix = createMatrix()
subMatrix = createMatrix()

bigSeq1, bigSeq2 = readInTrainingSet(infile1)
countDict, matchMatrix,dictionary, totalAA = compareSequences (matchMatrix,bigSeq1,bigSeq2,aminoAcidArray)
subMatrix = createSubMatrix(aminoAcidsT,matchMatrix,dictionary,totalAA,countDict)

print("This program uses a training set to create a substitution matrix. Please use sequences separated by new lines.")

stringSub = ''

print(aminoAcidL)
print(subMatrix)
outfile.write(">Training set substitution matrix \n")
outfile.write(aminoAcidL + '\n')
for i in range (1,21): #Traverse through subMatrix to create a string of the substitution index 
	for j in range(1,21):
		if j == 20:
			stringSub = stringSub + (str(subMatrix[i][j]) + '\n')
		else:
			stringSub = stringSub + (str(subMatrix[i][j]) + ',')

outfile.write(stringSub)
outfile.close()

