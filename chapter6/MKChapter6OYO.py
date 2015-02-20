'''
Name: Matthew Kachlik
Due: 2/10/15
Description: Implements Jukes-Cantor model 
'''
import math
import collections
from collections import defaultdict
'''
	Jukes
	measures evolutionary distance between two nucleotide sequences
'''
def Juke(seq1,seq2):

	difCtr = float(0)
	lenCtr = float(0)

	#count differences between sequences ignoring gaps
	for i in range(0, len(seq2)):
		if(seq1[i] != '-' and seq2[i] != '-'):
			lenCtr += 1
			if(seq1[i] != seq2[i]):
				difCtr += 1

	#Calculate and output results
	difCtr = difCtr/lenCtr
	jukes = (-3.0/4) * math.log(1-(4.0/3*difCtr))
	print (jukes)
	return jukes
def  Kimura(seq1,seq2,aDict):
	print (aDict)
	s = float(0)
	v = float(0)
	trans = float(0)
	sub = float(0)
	lenCtr = float(0)
	kimura = float(0)
	for i in range(0, len(seq2)):
		if(seq1[i] != '-' and seq2[i] != '-'):
				lenCtr += 1
				if(seq1[i] != seq2[i]):
					if(aDict[seq1[i]][seq2[i]] == 1):
						trans += 1

					elif(aDict[seq1[i]][seq2[i]] == 2):
						sub += 1


	s = float(trans/lenCtr)
	v = float(sub/lenCtr)
	kimura = 1.0/2 * math.log(1/(1-2*s-v)) + 1.0/4*math.log(1/(1-2*v))
	print (kimura)
def Tamura(seq1,seq2,aDict):
	s = float(0)
	v = float(0)
	trans = float(0)
	sub = float(0)
	lenCtr = float(0)
	tamura = float(0)
	gc1 = 0.0
	gc2 = 0.0
	c = float(0)
	for i in range(0, len(seq2)):
		
		if(seq1[i] != '-' and seq2[i] != '-'):
			lenCtr += 1
			if(seq1[i] == 'G' or seq1[i] == 'C'):
				gc1 += 1
			if(seq2[i] == 'G' or seq2[i] == 'C'):
				gc2 += 1
			if(seq1[i] != seq2[i]):
				if(aDict[seq1[i]][seq2[i]] == 1):
					trans += 1
				elif(aDict[seq1[i]][seq2[i]] == 2):
					sub += 1

	gc1 =  float(gc1/len(seq1))	
	gc2 =  float(gc2/len(seq2))	
	c = gc1 + gc2 - 2 * gc1 * gc2
	s = trans/lenCtr
	v = sub/lenCtr
	tamura = (c*-1.0)*math.log(1-s/c-v) - 1.0/2 *(1-c)*math.log(1-2*v)
	print (tamura)
def buildDirectionalString(matrix,N,M): 
    dString = ""
    currentRow = N
    currentCol = M
    gap = -1
    
    while(currentRow != 0 or currentCol != 0):
        if(currentRow == 0):
            dString =  dString + "H" # H is for Horizontal
            currentCol = currentCol-1
        elif(currentCol == 0):
            dString = dString + "V" # V is for Vertical
            currentRow = currentRow-1
        elif((matrix[currentRow][currentCol-1] + gap) == matrix[currentRow][currentCol]):
            dString = dString + "H"
            currentCol = currentCol-1
        elif((matrix[currentRow-1][currentCol] + gap) == matrix[currentRow][currentCol]):
            dString = dString + "V"
            currentRow = currentRow-1
        else:
            dString = dString + "D" # D is for Diagonal
            currentCol = currentCol-1
            currentRow = currentRow-1
    return dString

def semiGlobal(seq1,seq2):
	N = len(seq1) 
	M = len(seq2)

	match = int(raw_input("Please enter the match score "))
	misMatch = int(raw_input("Please enter the mismatch score "))

	matrix = [[0 for col in range(M+1)] for row in range(N+1)]
	
	gap = int(raw_input("Please enter the gap score "))
	for i in range (1,N+1):
	    matrix[i][0] = 0
	for i in range (1,M+1):
	    matrix[0][i] = 0
	for row in range(1,N+1):
	    for col in range (1,M+1):
	        if(seq1[row-1] == seq2[col-1]):
	            score1 = matrix[row-1][col-1] + match
	        else:
	            score1 = matrix[row-1][col-1] + misMatch
	        score2 = matrix[row][col-1] + gap
	        score3 = matrix[row-1][col] + gap
	        matrix[row][col] = max(score1,score2,score3)
	# Steps 2: Create strings that show direction of sequences:
	dString = buildDirectionalString(matrix,N,M)



	# Step 3: Build alignments using the directional strings:
	seq1Pos = N - 1 # position of last character in seq1
	seq2Pos = M - 1
	dirPos = 0

	#getting user data for scoring
	alignmentString = ""  
	tGapScore = 0
	matchScore = 0
	gapScore = 0
	#determing the alignment sequence based off the directional sequence
	while (dirPos < len(dString)):
	    if(dString[dirPos] == "D"):
	        if(seq1[seq1Pos] == seq2[seq2Pos]):
	            alignmentString = "|" + alignmentString
	            matchScore = matchScore + 1
	        else:
	            alignmentString = "*" + alignmentString
	        seq1Pos = seq1Pos - 1 
	        seq2Pos = seq2Pos - 1
	    elif(dString[dirPos] == "V"):
	        seq2List.insert(seq2Pos + 1, "-")
	        alignmentString = " " + alignmentString
	        seq1Pos = seq1Pos - 1
	        if seq2Pos == M-1 or seq2Pos == 0:
	            tGapScore += 1
	        else:
	            gapScore += 1
	    else:
	        seq1List.insert(seq1Pos + 1, "-")
	        alignmentString = " " + alignmentString
	        seq2Pos = seq2Pos - 1
	        if seq1Pos == N-1 or seq1Pos == 0:
	            tGapScore += 1
	        else:
	            gapScore += 1
	    dirPos = dirPos + 1
	return seq1,seq2
'''
	main
'''

dictBuild = "alignment.txt"
# Opens files
infile2 = open("giraffe.txt", 'r')
infile1 = open("whale.txt", 'r')
s1 = ""
# Reads the files into strings.
infile1.readline()  # bypass > header line
for line in infile1:
	line = line.replace('\n', '')
	s1 = s1 + line
s1 = s1.upper()
print (s1)

s2 = ""
infile2.readline()  # bypass > header line
for line in infile2:
    line = line.replace('\n', '')
    s2 = s2 + line
s2 = s2.upper()
infile2.close()
infile1.close()

alignDict = collections.defaultdict(dict)
alignDict = collections.defaultdict(lambda:defaultdict(int))

alignDict["G"]["A"] = 1
alignDict["A"]["G"] = 1
alignDict["C"]["T"] = 1
alignDict["T"]["C"] = 1

alignDict["A"]["C"] = 2
alignDict["C"]["A"] = 2
alignDict["G"]["T"] = 2
alignDict["T"]["G"] = 2
alignDict["A"]["T"] = 2
alignDict["T"]["A"] = 2
alignDict["C"]["G"] = 2
alignDict["G"]["C"] = 2
'''
infile3 = open(dictBuild,'r')
aminoArray = infile3.readline()
aminoArray = aminoArray.replace('\r','')
aminoArray = aminoArray.replace('\n','')
aminoList = list(aminoArray.split(","))
aminoAcidArray = [0 for row in range(len(aminoList))]

row = 0
for i in aminoList:
	if( row != len(aminoList) - 1):
		outfile.write(i+", ",)
	else:
		outfile.write(i)
	aminoAcidArray[row] = i
	row = row + 1

alignDict = collections.defaultdict(dict)
alignDict = collections.defaultdict(lambda:defaultdict(int))

for row in range(len(aminoList)):
	for col in range(len(aminoList)):
		if()
			alignDict[aminoAcidArray[row]][aminoAcidArray[col]] = float(1)
'''
#s1, s2 = semiGlobal(s1,s2)
while True:
	choice = input("Please enter which test you would like to run\n1) Jukes\n2) Kimura\n3)Tamura\nanything else program will quit\n")
	if(choice == "1"):
		Juke(s1,s2)
	elif(choice == "2"):
		Kimura(s1,s2,alignDict)
	elif(choice == "3"):
		Tamura(s1,s2,alignDict)
	else:
		break
