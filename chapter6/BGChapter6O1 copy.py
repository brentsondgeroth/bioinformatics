'''
Name: Brent Gaither
Description: This program uses the Jukes-Cantor, Tamura and Kimura models 
to estimate evolutionary distances.
Due Date: 2/17/15

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
    sequence = sequence.upper()
    return sequence

'''
compareSequences- reads in the sequences and compares to find matches and length of sequences.
'''

def compareSequences(alSeq1List,alSeq2List):
	lenCtr = 0
	difCtr = 0
	for char1, char2 in zip(alSeq1List, alSeq2List):
		if char1 != "-" and char2 != "-":
			lenCtr += 1
			if char1!= char2:
				difCtr += 1
	return lenCtr,difCtr

'''
buildMatrix- Builds the matrix for a semi global alignment
'''

def buildMatrix(gapScore,misMatch,match,N,M,sequence1,sequence2):

    #build matrix for semi global alignment
    matrix = [[0 for col in range(M+1)] for row in range(N+1)]
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

    return matrix
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
printMatrix- Uses the lists that make up the matrix and formats
to print in a table.
'''
def printMatrix(matrix):
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print ('\n'.join(table))
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

	alSeq1List = (alSeq1)
	alSeq2List = (alSeq2)
	return alSeq1List,alSeq2List

'''
findTransitions- Compares sequences to determine if the mutation was a 
transition or a transversion. 
'''

def findTransitions(alSeq1List,alSeq2List,lenCtr):
	key = ''
	matchCtr = 0
	gapCtr = 0
	transversion = 0
	transition = 0
	mutationType = {"GC":'transversion',"CG":'transversion', "AT":'transversion',"TA":'transversion',
	"AG":'transition',"GA":'transition',"CT":'transition',"TC":'transition',"GT":'transversion',"CA":'transversion',
	"AC":'transversion',"TG":'transversion'} #Creates a dictionary to determine relationship of mutation
	for char1, char2 in zip(alSeq1List, alSeq2List):
		if char1 != '-' and char2 != '-':
			if char1 == char2:
				matchCtr += 1
			else:
				key = char1 + char2
				if mutationType[key] == 'transition':
					transition +=1
				else:
					transversion +=1

	V = transversion/lenCtr
	S = transition / lenCtr

	return S,V

'''
findGC- Uses input sequences to find the amount of G and C in the sequences. 
'''

def findGC(alSeq1List, alSeq2List,lenCtr):

	gcSeq1 = 0
	gcSeq2 = 0
	seq1Len = 0
	seq2Len = 0
	for char1, char2 in zip(alSeq1List, alSeq2List):
		if char1 != '-': #IS THIS RIGHT OR LENCTR?
			seq1Len += 1
		if char2 != '-':
			seq2Len += 1
		if (char1 == "G" or char1 == "C") and char2 !="-": #Adds up the number of Gs and Cs in sequence
			gcSeq1 += 1
		if (char2 == "G" or char2 == "C") and char2 !="-":
			gcSeq2 += 1

	gcContent1 = gcSeq1 / lenCtr
	gcContent2 = gcSeq2 / lenCtr
	return gcContent1, gcContent2

'''
JukesCantor- Finds the evolutionary distance between sequences using 
Jukes-Cantor, a one parameter model.
'''

def findJukesCantor(difCtr,lenCtr):
	difCtr = difCtr / lenCtr #converts to fraction of substitutions 
	jukes = ((-3/4) * math.log(1-(4/3*difCtr)))
	return jukes

'''
Kimura- Finds the evolutionary distance between sequences using 
Kimura, a two parameter model.
'''

def findKimura(S,V,lenCtr):
	kimura = 0

	kimura = ((1/2) * math.log(1/(1-2*S-V))) + ((1/4) * math.log(1/(1-(2*V))))
	return kimura

'''
findTamura- Finds the evolutionary distance between sequences using 
Tamura, a two parameter model.
'''

def findTamura(S,V,lenCtr,alSeq1List,alSeq2List):

	gcContent1, gcContent2 = findGC(alSeq1List, alSeq2List,lenCtr)

	C = gcContent1 + gcContent2 - (2 * gcContent1 * gcContent2)

	tamura = -C * math.log(1 - (S/C) - V) - (((1/2)* (1-C)) * math.log(1-(2*V)))
	return tamura

'''
Main
'''

infile1 = open('hippo.txt', 'r')
infile2 = open('giraffe.txt', 'r')

sequence1 = readIn(infile1)
sequence2 = readIn(infile2)

N = len(sequence1)
M = len(sequence2)

gapScore = -1
misMatch = 0
match = 1
choice = ''
alignment = ''

print("This program aligns and calculates the evolutionary distance of two DNA sequences.")

alignment = input("Is the sequence aligned (Y/N) ")
alignment = alignment.upper()

if alignment == 'N': #aligns sequence if unaligned 
	matrix = buildMatrix(gapScore,misMatch,match,N,M,sequence1,sequence2)
	dString = buildDirectional(matrix,N,M,gapScore)
	alSeq1List,alSeq2List = analyzeMatrix(matrix,dString,N,M)

lenCtr, difCtr = compareSequences(sequence1,sequence2)
print(lenCtr)
print(difCtr)
S,V = findTransitions(sequence1,sequence2,lenCtr)

while(choice != '4'):
	choice = input("Please chose a distance calculation. Jukes(1) Kimura(2) Tamura(3) exit(4) ")
	if choice == '1':
		jukes = findJukesCantor(difCtr,lenCtr)
		print(jukes)
	elif choice == '2':
		kimura = findKimura(S,V,lenCtr)
		print(kimura)
	elif choice == '3':
		tamura = findTamura(S,V,lenCtr,sequence1,sequence2)
		print(tamura)



