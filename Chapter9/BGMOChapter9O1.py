'''
Brent Gaither Matt Obzera
This program finds possible genes in a sequence for eukaryotes and prokaryotes.

Due Date 3/2/2015
'''

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

def findProkaryote(sequence):
	genesFound = list()
	found = None
	prevPosStart = None
	posEndList = list()
	proCtr = 0
	j= 0
	threshholdShine = float(input("Input threshold for Shine Dalgarno site as percentage out of 100: "))/100
	threshholdPromoter = float(input("Input threshold for promoter as percentage out of 100: ")) /100
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

def findEukaryote(sequence):
	genesFound = list()
	found = None
	prevPosStart = None
	threshholdgcBox = float(input("Input threshold for gc box as percentage out of 100: ")) /100 #User input
	threshholdCatBox = float(input("Input threshold for cat box as percentage out of 100: ")) /100
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

#Main
stopCodons = {"TGA","TAA","TAG"}
infile = open("testSeq1.txt",'r')

sequence = readIn(infile)

organism = int(input("Kingdom of gene prokaryote(1) eukaryote(2): ")) #User inputs
orfSize = int(input("Input minimum size of ORF in nucleotides: "))

if orfSize % 3 != 0: #Ensure the open reading frame size stays in the same reading frame 
	orfSize = orfSize - orfSize % 3

#prokaryote finder 
if organism == 1:
	genesFound,found = findProkaryote(sequence)

#Eukaryokte finder
if organism == 2:
	genesFound,found = findEukaryote(sequence)
if found:
	for i in range(0,len(genesFound)):
		print(genesFound[i])
else:
	print("No genes found")


