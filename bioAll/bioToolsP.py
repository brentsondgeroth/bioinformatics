from __future__ import division
import collections
from collections import defaultdict
from collections import Mapping
import math
import sys
import re
import bioToolsM

'''
Description: Program that merges clusters using neighbor joining and outputs the final
product in Newick format with branch lengths.
'''
def phylogeneticTree():

	
	print("This program uses FASTA formatted sequences to create a phylogenetic tree in newick format using Jukes-Cantor as a distance calculator.")
	outfile = open("bgwmCh7skillsout.txt", 'w')# Read in data for nested hash structure
	
	userFile = raw_input ("Enter the file you wish to use: ")
	infile = open(userFile,'r')
	sequenceList,clusterNames = bioToolsM.readInList(infile)
	

	numClusters = len(clusterNames)# Set flag equal to the total number of clusters

	# Build nested hash structures
	originalDist = defaultdict(lambda:defaultdict(float))
	clusterDist = defaultdict(lambda:defaultdict(float))
	tempDist = defaultdict(lambda:defaultdict(float))

	#fill originalDist with distances calculated from Jukes-Cantor
	for i in range(0,int(numClusters)):
	    for j in range(0,int(numClusters)):
	        originalDist[clusterNames[i]][clusterNames[j]] = round((bioToolsM.findJukesCantor(sequenceList[i],sequenceList[j])),10)
	        clusterDist[clusterNames[i]][clusterNames[j]] = round((bioToolsM.findJukesCantor(sequenceList[i],sequenceList[j])),10)


	# Initialize a nested dictionary to build the Newick format output
	newick = defaultdict(lambda:defaultdict(float))

	while(numClusters > 2):
	    
	    transition,r = bioToolsM.calculateTransition(originalDist,clusterNames)
	    shortestD = bioToolsM.findSmallest(transition,clusterNames)

	    # Finds keys of smallest value
	    shortestK = bioToolsM.findKeys(transition, shortestD)
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
	    branchDis = bioToolsM.branchLength(originalDist,r,shortestI,shortestJ)
	    #Add in the distances for the branches with branchDis
	    myValue = '(' + x +":" + str(branchDis[shortestI]) + ',' + y +":" + str(branchDis[shortestJ]) + ')'
	  
	    newick[myKey][myValue] = 0
	    
	    # merge clusters I and J
	    newCluster = shortestI + shortestJ

	    for iName  in clusterNames: #Set the original distances to the new distances 
	        for jName  in clusterNames:
	            originalDist[iName][jName] = clusterDist[iName][jName]

	    #Calculate neighbor joining distances
	    clusterDist = bioToolsM.neighborJoining(originalDist,newCluster, clusterNames,shortestI,shortestJ,clusterDist)

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
Description:This program finds possible genes in a sequence for eukaryotes and prokaryotes.
'''
def geneProbablity ():
	stopCodons = {"TGA","TAA","TAG"}
	
	userFile = raw_input ("Enter the file you wish to use: ")
	infile = open(userFile,'r')

	sequence = bioToolsM.readIn(infile)

	organism = int(raw_input("Kingdom of gene prokaryote(1) eukaryote(2): ")) #User inputs
	orfSize = int(raw_input("Input minimum size of ORF in nucleotides: "))

	if orfSize % 3 != 0: #Ensure the open reading frame size stays in the same reading frame 
		orfSize = orfSize - orfSize % 3

	#prokaryote finder 
	if organism == 1:
		genesFound,found = bioToolsM.findProkaryote(sequence,orfSize,stopCodons)

	#Eukaryokte finder
	if organism == 2:
		genesFound,found = bioToolsM.findEukaryote(sequence,orfSize,stopCodons)
	if found:
		for i in range(0,len(genesFound)):
			print(genesFound[i])
	else:
		print("No genes found")
'''
Description: This program uses a Hidden Markov model to calculate the probability of a 
eukaryotic sequence. 
'''
def markovModel ():
	print("Find the probability of a eukaryotic gene from a FASTA file.")
	userFile = raw_input ("Enter the file you wish to use: ")
	infile = open(userFile,'r')

	sequence = bioToolsM.readIn(infile)

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
	bioToolsM.markov(states,'B',length,ctr,'X','',pathList)

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
	print("Align DNA sequences from two FASTA files")
	userFile1 = raw_input ("Enter the first file you wish to use: ")
	userFile2 = raw_input ("Enter the second file you wish to use: ")
	infile1 = open(userFile,'r')
	infile2 = open(userFile,'r')

	outfile = open('needlemanOut.txt', 'w')
	s1 = bioToolsM.readIn(infile1)
	s2 = bioToolsM.readIn(infile2)

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

	bioToolsM.buildMatrix(matrix,alignment,gap,misMatch,match,N,M,s1,s2)
	#step 2 create directional string
	dString = ''
	maxValue = ''
	if(alignment == 'G' or alignment == 'S'):
	    dString = bioToolsM.buildDirectional(matrix, N, M, gap,alignment)
	else:   
	    maxValue = max([max(row) for row in matrix]) #Find the best alignment in matrix
	    for rowIdx, rowValue in enumerate(matrix):
	        if max(rowValue) == maxValue:
	            maxRowIdx = rowIdx
	            maxColIdx = rowValue.index(max(rowValue))
	            N = maxRowIdx
	            M = maxColIdx
	            dString = bioToolsM.buildDirectional(matrix, N, M, gap,alignment)
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
		bioToolsM.printMatrix(matrix,outfile)
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

	userFile = raw_input ("Enter the file for the training set you wish to use: ")
	infile = open(userFile,'r')
	outfile = open('substitutionMatrix.txt', 'w')

	aminoAcids  = ("0,A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V") # amino acids used to create matrix
	aminoAcidL = ("A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V")#amino acid list used for out file
	aminoAcidArray = aminoAcids.split(',')
	aminoAcidsT = ('A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'); #tuple of amino acids used for dict

	matchMatrix = bioToolsM.createMatrix(aminoAcidArray)
	subMatrix = bioToolsM.createMatrix(aminoAcidArray)

	bigSeq1, bigSeq2 = bioToolsM.readInTrainingSet(infile1)
	countDict, matchMatrix,dictionary, totalAA = bioToolsM.compareSequences(matchMatrix,bigSeq1,bigSeq2,aminoAcidArray)
	subMatrix = bioToolsM.createSubMatrix(aminoAcidsT,matchMatrix,dictionary,totalAA,countDict,subMatrix)

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