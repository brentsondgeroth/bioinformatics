'''
Brent Gaither 
Program that combines multiple bioinformatics programs into one large program

'''

from __future__ import division
import collections
from collections import defaultdict
from collections import Mapping
import math
import sys
import re

'''
readIn- Reads in a FASTA formatted file 
'''

def readIn (infile):
	infile.readline()#Read first line of Fasta formatted filed
	sequence = ""
	for line in infile:#Read in remaining data in Fasta file
		line.replace("\n","")
		line.replace("\r","")
		sequence = sequence + line
	sequence = sequence.upper()
	return sequence

'''
readInList- Reads in a fasta formatted file into a list 
'''

def readInList(infile):
    data=''
    sequenceNames=[]
    sequenceList=[]
    line =infile.readline()

    sequenceNames.append(line[1])
    for line in infile:
        if line.startswith('>'):
            sequenceNames.append(line[1])
            sequenceList.append(data)
            data=''
        else:
            data= data + line.upper().strip()
    sequenceList.append(data)
    return sequenceList,sequenceNames
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
findPattern- Takes in a sequence and finds if the pattern is in it. Uses a threshold to
determine if the pattern is close. 
'''

def findPattern(pattern,searchText,startLoc,stopLoc,increment,threshhold):
	textLen = len(searchText)
	patternLen = len(pattern)
	if stopLoc > textLen: #Ensures stop location is within the size of the sequence
		stopLoc = textLen
	if startLoc < 0: #Ensures start is within size of sequence
		startLoc = 0
	if startLoc < stopLoc:	
		for i in range(startLoc,stopLoc,increment):
			ctr= 0
			j = i
			for k in range(0,patternLen):
				if searchText[j] == pattern[k]:
					ctr +=1 #Counts matches of a pattern
				j+=1
			if (float(ctr)/float(patternLen)) >= threshhold:
				return i
	return -1

'''
findProkaryote- Method to find the possible genes in a prokaryotic sequence.
'''

def findProkaryote(sequence,orfSize,stopCodons):
	genesFound = list()
	found = None
	prevPosStart = None
	posEndList = list()
	proCtr = 0
	j= 0
	threshholdShine = float(raw_input("Input threshold for Shine Dalgarno site as percentage out of 100: "))/100
	threshholdPromoter = float(raw_input("Input threshold for promoter as percentage out of 100: ")) /100
	while len(sequence) >= 80 + proCtr:#Find all sequences in the entire sequence
		pattern = "ATG"#Find start codon
		posStart = findPattern(pattern,sequence,80 + proCtr,len(sequence) - orfSize + 3,1,1)
		if posStart != -1:
			posEnd = -1
			posEndList = []
			for stop in stopCodons:#Finds all stop codons
				tempEnd = (findPattern(stop,sequence,posStart + 1 + orfSize,len(sequence)-3,3,1))#Traverse to find stop codon in same reading frame
				if tempEnd > posStart:
					posEndList.append(tempEnd)
			if posEndList != []:
				posEnd = min([x for x in posEndList if x !=-1])#finds the first stop codon in the reading frame
			if posEnd != -1:			
				pattern = "AGGAGG"
				shineDalgarno = findPattern(pattern,sequence,posStart - 15,posStart-3,1,threshholdShine)
				if shineDalgarno != -1:
					pattern = "TATAAT" #Check for -10 
					negTen = findPattern(pattern,sequence,posStart-500,posStart,1,threshholdPromoter)
					if negTen != -1:
						j +=1
						pattern = "TTGACA" #Check for -35 promoter
						negThirtyFive = findPattern(pattern,sequence,negTen-30,negTen,1,threshholdPromoter)
						if negThirtyFive != -1 and posStart != prevPosStart:
							found = True
							prevPosStart = posStart#Stops from adding the same gene twice
							genesFound.append("Possible gene found at: " + str(posStart) + " ending at " + str(posEnd))
		proCtr+=1
	return genesFound, found

'''
findEukaryote- Method to find the possible genes in a eukaryotic sequence.
'''

def findEukaryote(sequence,orfSize,stopCodons):
	genesFound = list()
	found = None
	prevPosStart = None
	threshholdgcBox = float(raw_input("Input threshold for gc box as percentage out of 100: ")) /100 #User input
	threshholdCatBox = float(raw_input("Input threshold for cat box as percentage out of 100: ")) /100
	eukCtr = 0
	while len(sequence) >= 80 + eukCtr:#Find all sequences in the entire sequence
		pattern = "ATG"
		posStart = findPattern(pattern,sequence,80 + eukCtr,len(sequence) - orfSize + 3,1,1)
		if posStart != -1:
			tataSegment = sequence[posStart-min(500,posStart):posStart]#Takes only the section of the DNA needed for TATA box
			tataBox = re.search(r'TATA[AT]A[AT]', tataSegment) 
			if tataBox == None:
				tataStart = -1
			else:
				tataStart = tataBox.start()
			if tataStart != -1:
				inrSegment = sequence[posStart-min(500,posStart):posStart]
				inr = re.search(r'[TC]{2}CA[AG]{2}',inrSegment)#regular expression
				if inr != None:#Ensure the regular expression returned a position
					inrStart = inr.start()
				else:
					inrStart = -1
				if inrStart != -1:
					pattern = "GGGCGG"
					gcBox = findPattern(pattern,sequence,tataStart-min(500, tataStart),tataStart,1,threshholdgcBox) 
					if gcBox != -1:
						pattern = "CAAT"
						catBox = findPattern(pattern,sequence,tataStart-500,tataStart,1,threshholdCatBox)
						if catBox != -1:
							kozakSegment = sequence[posStart-11:posStart+4]
							kozak = re.search(r'[GC]{3}[AG][ATGC]{2}ATGG',kozakSegment)
							if kozak != None:
								kozakStart = kozak.start()
							else:
								kozakStart = -1
							if kozakStart != -1 and posStart != prevPosStart:
								prevPosStart = posStart#stops from adding the same gene twice
								found = True
								genesFound.append("Possible gene found at: " + str(posStart))
		eukCtr +=1
	return genesFound, found



'''
markov- recursive method to create all the possible paths for 
the sequence
'''

def markov (states,currentState, length, ctr,end,path,pathList):

	if ctr == length and currentState == end:#End with base case of path with length "length" and ending with "end"
		pathList.append(path)
	elif currentState == end: #Return if go to end state but not with expected length
		return
	elif ctr<length: #Continue to recursively call function until path is made
		for state in states[currentState]:#Go through each possible state to create each unique path		
			markov(states,state, length, ctr+1,end,path+state,pathList)
	else:
		return

'''
findJukesCantor- Finds the evolutionary distance between sequences using 
Jukes-Cantor, a one parameter model.
'''

def findJukesCantor(alSeq1,alSeq2):
    lenCtr = 0
    difCtr = 0
    if alSeq1 == alSeq2:
        return 0
    for char1, char2 in zip(alSeq1, alSeq2): #traverses through sequences
        if char1 != "-" and char2 != "-":
            lenCtr += 1
            if char1!= char2:
                difCtr += 1
    difCtr = difCtr / lenCtr #converts to fraction of substitutions 
    jukes = ((-3.0/4) * math.log(1-(4.0/3*difCtr)))
    return jukes



'''
findKeys- Used the smallest value and the reverse copy of the nested dictionary
to find the value and return it's parent and child keys
'''

def findKeys(nested,value):
    for k,v in nested.items():
        if isinstance(v,dict):
            p = findKeys(v,value)
            if p:
                return [k] + p
        elif v == value:
            return [k]

'''
findSmallest- Finds and returns the shortest distance between sequences
listed in the matrix
'''

def findSmallest(transition,clusterNames):

    smallest = 10000
    for iName  in clusterNames:
        for jName  in clusterNames:
            if iName != jName:
                if smallest > float(transition[iName][jName]) and float(transition[iName][jName]) != 0:
                    smallest = float(transition[iName][jName])
    return smallest
       
'''
singleLinkage- Calculate distances between new cluster and all
other clusters using single linkage
'''

def singleLinkage(clusterDist, newCluster, originalDist, clusterNames):
    for cluster in clusterNames:
        smallestD = max(clusterDist)
        for c1 in cluster:
            for c2 in newCluster:
                if originalDist[c1][c2] < smallestD:
                    smallestD = originalDist[c1][c2]
        clusterDist[newCluster][cluster] = smallestD
        clusterDist[cluster][newCluster] = smallestD
    return clusterDist

'''
calculateTransition- Calculates the transition matrix for a neighbor joining phylogeny.
'''

def calculateTransition(originalDist,clusterNames):
    r = collections.defaultdict(float)
    transition = defaultdict(lambda:defaultdict(float))
    dix = 0
    i =0
    j = 0
    for iName  in clusterNames:
        for jName  in clusterNames:
            dix = dix + int(originalDist[iName][jName]) #calculate distance between i and j 
            j += 1
        r[iName] = dix / (len(clusterNames) - 2) #Calculates the r values for the transition matrix
        i +=1
        dix  = 0
    k =0
    l = 0
    for iName  in clusterNames:
        for jName  in clusterNames:
            transition[iName][jName] = originalDist[iName][jName]  - r[iName] - r[jName] #Fill transition matrix with transition distances 
            k +=1
            l+=1
    return transition, r

'''
neighborJoining- Calculates the neighbor joining distances. Does not assume constant rate of evolution/change 
'''

def neighborJoining(originalDist,newCluster, clusterNames,shortestI,shortestJ,clusterDist):
    for iName  in clusterNames:
        for jName  in clusterNames:
            if iName == jName:
                clusterDist[iName][jName] = 0
            else:
                clusterDist[newCluster][jName] = (float(originalDist[shortestI][iName])) + float(originalDist[shortestJ][jName]) - (float(originalDist[newCluster[0]][newCluster[1]]))/2
    return clusterDist
    
'''
branchLength- Determines distance from two clusters and their most recent ancestor.
'''

def branchLength(originalDist,r,shortestI,shortestJ):
    branchDis = collections.defaultdict(float)
    branchDis[shortestI] = (originalDist[shortestI][shortestJ] + r[shortestI] - r[shortestJ])/2
    branchDis[shortestJ]= (originalDist[shortestI][shortestJ] + r[shortestJ] - r[shortestI])/2
    return branchDis

'''
printMatrix- Uses the lists that make up the matrix and formats
to print in a table.
'''
def printMatrix(matrix,outfile):
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print ('\n'.join(table))
    outfile.write('\n'.join(table))
'''
buildMatrix- Builds a string that provides direction for tracing
a path back through the matrix to build the best match.
'''
def buildMatrix(matrix,alignment,gap,misMatch,match,N,M,s1,s2):

    #Step 1 build matrix
    matrix [0][0] = 0
    for i in range(1,N+1): 
        if alignment == 'G':
            matrix[i][0] = matrix[i-1][0]+int(gap) #Sets col to penalize for initial or terminal gaps
        else:
            matrix [i][0] = matrix[i-1][0] #No penalty for local or semi local initial or terminal gaps
        for j in range(1,M+1):
            if alignment == 'G':
                matrix [0][j] = matrix[0][j - 1] + int(gap) 
            else:
                matrix [0][j] = matrix[0][j-1] 
            if (s1[i-1] == s2[j - 1]):
                score1 = matrix[i - 1][j - 1] + int(match)               
            else:
                score1 = matrix[i - 1][j - 1] + int(misMatch)
            if alignment == 'G':
                score2 = matrix[i][j - 1] + int(gap)
                score3 = matrix[i - 1][j] + int(gap)
            else:
                score2 = matrix[i][j - 1]
                score3 = matrix[i-1][j]
            if (alignment == 'L'):
                matrix[i][j] = max(score1, score2, score3, 0) #Does not allow for negative numbers in local alignment
            else:
                matrix[i][j] = max(score1, score2, score3) 

    return matrix
'''
buildDirectional- Creates the directional string for following the matrix to 
the best solution
'''
def buildDirectional(matrix, N, M, gapScore,alignment):
    dString = ''
    currentRow = N
    currentCol = M
    while(currentRow != 0 or currentCol != 0):

        if(alignment == 'L' and matrix[currentRow][currentCol] == 0): #returns dstring if local alignment is complete
            return dString
        if(currentRow == 0):
            dString = dString + ('H')
            currentCol = currentCol - 1
        elif(currentCol == 0):
            dString = dString + ('V')
            currentRow = currentRow - 1
        elif(matrix[currentRow][currentCol] == matrix[currentRow-1][currentCol] and alignment == 'S' and currentCol == M):
            dString = dString + ('V') #Moves semi local alignment away from terminal gaps
            currentRow = currentRow-1
        elif(matrix[currentRow][currentCol] == matrix[currentRow][currentCol-1] and alignment == 'S' and currentRow == N):
            dString = dString + ('H')
            currentCol = currentCol -1
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

def createMatrix(aminoAcidArray):
	matrix = [[1 for col in range(21)] for row in range(21)]#Fills matrix with 1's to initialize and uses a pseudo count
	for i in range(0,21):
		matrix[i][0] = aminoAcidArray[i]
		for j in range(0,21):
			matrix[0][j] = aminoAcidArray[j]
	return matrix
'''
createSubMatrix- Creates and computes calculations for the substitution matrix. 
'''
def createSubMatrix(aminoAcidsT,matchMatrix,dictionary,totalAA,countDict,subMatrix):
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
Program that merges clusters using neighbor joining and outputs the final
product in Newick format with branch lengths.
'''
def phylogeneticTree():

	# Read in data for nested hash structure
	outfile = open("bgwmCh7skillsout.txt", 'w')
	infile = open("file3.txt",'r')
	sequenceList,clusterNames = readInList(infile)
	print("This program uses FASTA formatted sequences to create a phylogenetic tree in newick format using Jukes-Cantor as a distance calculator.")

	numClusters = len(clusterNames)# Set flag equal to the total number of clusters

	# Build nested hash structures
	originalDist = defaultdict(lambda:defaultdict(float))
	clusterDist = defaultdict(lambda:defaultdict(float))
	tempDist = defaultdict(lambda:defaultdict(float))

	#fill originalDist with distances calculated from Jukes-Cantor
	for i in range(0,int(numClusters)):
	    for j in range(0,int(numClusters)):
	        originalDist[clusterNames[i]][clusterNames[j]] = round((findJukesCantor(sequenceList[i],sequenceList[j])),10)
	        clusterDist[clusterNames[i]][clusterNames[j]] = round((findJukesCantor(sequenceList[i],sequenceList[j])),10)


	# Initialize a nested dictionary to build the Newick format output
	newick = defaultdict(lambda:defaultdict(float))

	while(numClusters > 2):
	    
	    transition,r = calculateTransition(originalDist,clusterNames)
	    shortestD = findSmallest(transition,clusterNames)

	    # Finds keys of smallest value
	    shortestK = findKeys(transition, shortestD)
	    shortestK = sorted(shortestK)
	    shortestI = shortestK[0]
	    shortestJ = shortestK[1]

	    '''
	    For each cluster pair of a given shortest distance,
	    merge them in Newick format and place them in newick as child keys.
	    If they were already in the dictionary individually, delete them.
	    '''
	    if shortestI in newick:
	        x = newick[shortestI].iterkeys().next()
	        del newick[shortestI]
	    else:
	        x = shortestI
	    if shortestJ in newick:
	        y = newick[shortestJ].iterkeys().next()
	        del newick[shortestJ]
	    else:
	        y = shortestJ

	    '''
	    The merged clusters become the parent key, and the merged newick clusters
	    become the child key.
	    '''
	    myKey = shortestI + shortestJ
	    branchDis = branchLength(originalDist,r,shortestI,shortestJ)
	    #Add in the distances for the branches with branchDis
	    myValue = '(' + x +":" + str(branchDis[shortestI]) + ',' + y +":" + str(branchDis[shortestJ]) + ')'
	  
	    newick[myKey][myValue] = 0
	    
	    # merge clusters I and J
	    newCluster = shortestI + shortestJ

	    for iName  in clusterNames: #Set the original distances to the new distances 
	        for jName  in clusterNames:
	            originalDist[iName][jName] = clusterDist[iName][jName]

	    #Calculate neighbor joining distances
	    clusterDist = neighborJoining(originalDist,newCluster, clusterNames,shortestI,shortestJ,clusterDist)

	    #remove the clusters that have been combined 
	    clusterNames.remove(shortestI)
	    clusterNames.remove(shortestJ)

	    #remove the child key and values
	    for k,v in clusterDist.items():
	        for k2,v2 in v.items():
	            if k2 == shortestI:
	                del clusterDist[k][k2]
	            if k2 == shortestJ:
	                del clusterDist[k][k2]

	    clusterNames.append(newCluster)

	    '''
	    Finally, remove parent keys and remaining nested keys of the
	    individual clusters after they have merged
	    '''
	    for k,v in clusterDist.items():
	        if k == shortestI:
	            del clusterDist[k]
	        if k == shortestJ:
	            del clusterDist[k]

	    print ("merging clusters " + shortestI + " and " + shortestJ)
	    outfile.write("merging clusters " + shortestI + " and " + shortestJ + '\n')
	    numClusters -= 1

	# Put child keys from Newick in list, then join list elements as a string
	nFormat = []

	for k,v in newick.items():
	    for k2,v2 in v.items():
	        nFormat.append(k2) 

	nFormat = ','.join(nFormat)
	print("merging clusters " + clusterNames[0] + " and " +clusterNames[1])
	# Surround the newick string with a final set of parentheses
	print ("((" + nFormat + ":" + str(branchDis[clusterNames[0]]) + ")" )
	outfile.write(','.join(clusterNames))
	outfile.write("\n((" + nFormat + ":" + str(branchDis[clusterNames[0]]) + ")")
	outfile.close()

'''
This program finds possible genes in a sequence for eukaryotes and prokaryotes.
'''

def geneProbablity ():
	stopCodons = {"TGA","TAA","TAG"}
	infile = open("testSeq1.txt",'r')

	sequence = readIn(infile)

	organism = int(raw_input("Kingdom of gene prokaryote(1) eukaryote(2): ")) #User inputs
	orfSize = int(raw_input("Input minimum size of ORF in nucleotides: "))

	if orfSize % 3 != 0: #Ensure the open reading frame size stays in the same reading frame 
		orfSize = orfSize - orfSize % 3

	#prokaryote finder 
	if organism == 1:
		genesFound,found = findProkaryote(sequence,orfSize,stopCodons)

	#Eukaryokte finder
	if organism == 2:
		genesFound,found = findEukaryote(sequence,orfSize,stopCodons)
	if found:
		for i in range(0,len(genesFound)):
			print(genesFound[i])
	else:
		print("No genes found")
'''
This program uses a Hidden Markov model to calculate the probability of a 
eukaryotic sequence. 
'''
def markovModel ():
	infile = open("file1.txt",'r')

	sequence = readIn(infile)

	#Create nested hash tables
	transition = defaultdict(lambda:defaultdict(float))
	emission = defaultdict(lambda:defaultdict(float))
	states = defaultdict(lambda:defaultdict())

	'''
	B = begin, W = start of sequence, A = alpha Exon, E = internal exon, S,P splice donor sites,
	C,R splice acceptor site, I = internal intron FGH = first ATG, XYZ = terminal nucleotides of terminal exon
	'''
	states = {'B':['W'],'W':['W','F'],'A':['A','S'],'F':['G'],'G':['H'],'H':['A'],'E':['E','X','S'],'S':['P'],'P':['I'],'I':['I','C'],'C':['R'],'R':['E'],'X':['X']}

	#Transition values
	transition['W']['W'] = .9
	transition['W']['F'] = .1
	transition['F']['G'] = 1
	transition['G']['H'] = 1
	transition['H']['A'] = 1

	transition['A']['A'] = .9
	transition['A']['S'] = .1

	transition['E']['E'] = .8
	transition['E']['S'] = .1
	transition['E']['X'] = .1
	transition['S']['P'] = 1
	transition['P']['I'] = 1

	transition['I']['I'] = .9
	transition['I']['C'] = .1
	transition['C']['R'] = 1
	transition['R']['E'] = 1

	#emission values 
	emission['C']['A'] = .9998
	emission['C']['T'] = .000067
	emission['C']['G'] = .000067
	emission['C']['C'] = .000067
	emission['R']['A'] = .0005
	emission['R']['T'] = .0001
	emission['R']['G'] = .9993
	emission['R']['C'] = .0001

	emission['S']['A'] = .0005
	emission['S']['T'] = .0001
	emission['S']['G'] = .9993
	emission['S']['C'] = .0001
	emission['P']['A'] = .0001
	emission['P']['T'] = .0069
	emission['P']['G'] = .0001
	emission['P']['C'] = .9929

	emission['F']['A'] = 1
	emission['G']['T'] = 1
	emission['H']['G'] = 1

	emission['E']['A'] = .2
	emission['E']['T'] = .3
	emission['E']['G'] = .3
	emission['E']['C'] = .2

	emission['I']['A'] = .27
	emission['I']['T'] = .3
	emission['I']['G'] = .23
	emission['I']['C'] = .2


	emission['A']['A'] = .2
	emission['A']['T'] = .3
	emission['A']['G'] = .3
	emission['A']['C'] = .2

	emission['W']['A'] = .2
	emission['W']['T'] = .3
	emission['W']['G'] = .3
	emission['W']['C'] = .2

	emission['X']['A'] = 1
	emission['X']['T'] = 1
	emission['X']['G'] = 1
	emission['X']['C'] = 1


	ctr = 0
	pathList = list()
	total = list()
	length = len(sequence)
	markov(states,'B',length,ctr,'X','',pathList)

	p = [float(1)] *(len(pathList))#Create list to hold probabilities 
	j = 0
	for path in pathList:
		previousState = 'W'#Set first previous state to sequence start
		i = 0
		for state in path:#multiply emission and transition values for total probability of a path
			total.append(float(transition[previousState][state]) * float(emission[state][sequence[i]]))
			previousState = state
			i+=1
		for num in total:
			p[j] *= float(num)
		j+=1
		del total[:]

	k = 0
	for num in p:#Calculate log of all probabilities
		if num!=0:
			p[k] = math.log(p[k])
		else:
			p[k] = 0
		k+=1
	'''
	for i in range(0,len(pathList)):#Prints all possible paths
		print(pathList[i])
	'''

	print("W = start of sequence, A = alpha Exon, E = internal exon, S,P splice donor sites," +
		"C,R splice acceptor site, I = internal intron FGH = first ATG, XYZ = terminal nucleotides of terminal exon")
	print("Top sequence is the path list for the most likely eukaryotic start sequence (bottom)")

	line = 0
	nucleotideCount = 1
	charPerLine = 50
	startSite = max([x for x in p if x !=0])

	# Prints the alignment with a limit of 50 characters per line.
	for line in range(line,len(sequence),charPerLine):
	    print(str(nucleotideCount) + " " + pathList[(p.index(startSite))][line:line + charPerLine])
	    print(str(nucleotideCount) + " " + sequence[line:line+charPerLine])
	    nucleotideCount = nucleotideCount + charPerLine
'''
Description: This program uses the needleman wunch algorithm to align 
sequences globally or semi globally or locally 
'''
def needleman ():

	# Opens files
	infile1 = open('file1.txt', 'r')
	infile2 = open('file2.txt', 'r')
	outfile = open('ch3sk3out.txt', 'w')
	s1 = readIn(infile1)
	s2 = readIn(infile2)

	temp = ""
	if (len(s1) > len(s2)): #ensure long sequence is s2
	    temp = s1
	    s1 = s2
	    s2 = temp

	#Initialize matrix and prompt use for score preferences
	N = len(s1)
	M = len(s2)
	connectors = ''
	matrix = ''
	matrix = list(matrix)
	matrix = [[0 for col in range(M+1)] for row in range(N+1)]
	print("This program aligns DNA sequences with global and semi global options.")
	gap = raw_input("Please enter gap score: ")
	misMatch = raw_input("Please enter mismatch score: ")
	match = raw_input("Please enter match score: ")
	alignment = raw_input("[G]lobal or [S]emi global alignment or [L]ocal: ")
	alignment = alignment.upper()

	if alignment == 'G':
		alignmentType = "Globally"
	elif alignment == 'S':
		alignmentType = "Semi Globally"
	else:
		alignmentType = "Locally"

	buildMatrix(matrix,alignment,gap,misMatch,match,N,M,s1,s2)
	#step 2 create directional string
	dString = ''
	maxValue = ''
	if(alignment == 'G' or alignment == 'S'):
	    dString = buildDirectional(matrix, N, M, gap,alignment)
	else:   
	    maxValue = max([max(row) for row in matrix]) #Find the best alignment in matrix
	    for rowIdx, rowValue in enumerate(matrix):
	        if max(rowValue) == maxValue:
	            maxRowIdx = rowIdx
	            maxColIdx = rowValue.index(max(rowValue))
	            N = maxRowIdx
	            M = maxColIdx
	            dString = buildDirectional(matrix, N, M, gap,alignment)
	#step 3 Build alignment using directional strings
	seq1Pos = N-1
	seq2Pos = M-1
	dirPos = 0
	alSeq1 = ''
	alSeq2 = ''
	alignmentScore = 0
	matrixOut = ''

	while(dirPos < len(dString)):
	        
	    if(dString[dirPos] == "D"): #Align sequence for match
	        alSeq1 = s1[seq1Pos] + alSeq1
	        alSeq2 = s2[seq2Pos] + alSeq2
	        seq1Pos = seq1Pos - 1
	        seq2Pos = seq2Pos - 1
	    elif(dString[dirPos] == "V"):#Place a space for a gap
	        alSeq1 = s1[seq1Pos] + alSeq1
	        alSeq2 = '-' + alSeq2
	        seq1Pos = seq1Pos - 1
	    else:
	        alSeq1 = '-' + alSeq1 #Place a space for a gap
	        alSeq2 = s2[seq2Pos] + alSeq2
	        seq2Pos = seq2Pos -1

	    dirPos = dirPos + 1

	alSeq1List = list(alSeq1)
	alSeq2List = list(alSeq2)

	for char1, char2 in zip(alSeq1List, alSeq2List):
	    if char1 == char2:
	        connectors = connectors + "|" #Create list for relationship between sequences
	        alignmentScore = alignmentScore + 1

	    else:
	        if char1 == "-" or char2 == "-":
	            connectors = connectors + " "
	        else:
	            connectors = connectors + "."

	line = 0
	nucleotideCount = 1
	charPerLine = 50

	# Prints the alignment with a limit of 50 characters per line.
	for line in range(line,len(alSeq1),charPerLine):
	    print(str(nucleotideCount) + " " + alSeq1[line:line + charPerLine])
	    outfile.write(str(nucleotideCount) + " " + alSeq1[line:line + charPerLine]+ '\n')
	    print(len(str(nucleotideCount)) * " " + " " + connectors[line:line + charPerLine])
	    outfile.write(len(str(nucleotideCount)) * " " + " " + connectors[line:line + charPerLine]+ '\n')
	    print(str(nucleotideCount) + " " + alSeq2[line:line+charPerLine])
	    outfile.write(str(nucleotideCount) + " " + alSeq2[line:line+charPerLine])
	    nucleotideCount = nucleotideCount + charPerLine
	print("Alinged " + alignmentType)
	print("Match percentage is: " + str(alignmentScore/float(len(alSeq1))*100))
	outfile.write("\nMatch percentage is: " + str(alignmentScore/float(len(alSeq1))*100) + '\n')
	print("The alignment score is: " + str(matrix[N][M]))
	outfile.write("The alignment score is: " + str(matrix[N][M])+ '\n' )
	matrixOut = raw_input("Would you like to print the matrix? (1 yes, 2 no) ")
	if (matrixOut == "1"):
		printMatrix(matrix,outfile)
	dStringOut = raw_input("Would you like to print the path string? (1 yes, 2 no) ")
	if (dStringOut == "1"):
	    print (dString)
	    outfile.write('\n' + dString)
'''
Description: This program uses a training set to create a substitution matrix. Use sequences 
separated by new lines to differentiate sequences. Do not include headers for sequences. 
Place sequences one after the other to be read in as pairs. 
'''
def substitutionMatrix():

	infile1 = open('ecoliTrainSet.txt', 'r')
	outfile = open('ch5O1out.txt', 'w')

	aminoAcids  = ("0,A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V") # amino acids used to create matrix
	aminoAcidL = ("A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V")#amino acid list used for out file
	aminoAcidArray = aminoAcids.split(',')
	aminoAcidsT = ('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'); #tuple of amino acids used for dict

	matchMatrix = createMatrix(aminoAcidArray)
	subMatrix = createMatrix(aminoAcidArray)

	bigSeq1, bigSeq2 = readInTrainingSet(infile1)
	countDict, matchMatrix,dictionary, totalAA = compareSequences (matchMatrix,bigSeq1,bigSeq2,aminoAcidArray)
	subMatrix = createSubMatrix(aminoAcidsT,matchMatrix,dictionary,totalAA,countDict,subMatrix)

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

