'''
Brent Gaither
	
This program creates fragments of DNA from a Fasta formatted file. The fragments
are then overlapped to recreate the original fragment using the greedy approach for 
the traveling sales person problem.
Due Date 2/24/15
'''
from random import randint
from collections import defaultdict
import operator


'''
findBiggest- Find the biggest value in a nested dictionary
'''

def findBiggest(dictionary,frags):
    biggest = 0
    for iName  in range(0,len(frags)-1):
        for jName in range(0,len(frags)-1):
                if biggest < int(dictionary[iName][jName]):
                    biggest = int(dictionary[iName][jName])
    return biggest

'''
findKeys- Find the value and return it's parent and child keys
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
coverageMet- Determines if coverage for the sequence is large enough
'''

def coverageMet(coverage,fold):
	i = 0
	met = True #Assume coverage met
	while(i<len(coverage) and met):
		if coverage[i] < fold:
			met = False
		i += 1
	return met

'''
createDict- Iterates through the frags calling the findOverlap method 
to help create the dictionary of aligned sequences.
'''

def createDict(frags):

	#Creates the dictionary of fragments and its overlaps with each other the fragments
	fragDict = defaultdict(lambda:defaultdict(int))
	numFrags = len(frags)
	for i in range(0,numFrags-1):
		for j in range(i+1,numFrags):
			frag1 = frags[i]
			frag2 = frags[j]
			contig, overlap, checker = findOverlap(frag1,frag2)
			if checker: #Determines the position of the fragment in the dictionary.
				fragDict[frags.index(frag1)][frags.index(frag2)] = overlap
			else:
				fragDict[frags.index(frag2)][frags.index(frag1)] = overlap
	return fragDict

'''
findOverlap- Finds the overlap between two sequences and returns the combined sequence(contig)
and the overlap.
'''

def findOverlap(frag1,frag2):

	f1Len = len(frag1)
	f2Len = len(frag2)
	minLen = min(f1Len,f2Len)
	overlap = 0
	k = minLen -1 
	contig = ""
	checker = None
	while k >= 1 and overlap == 0: #check for overlap of sequences 
		if frag1[f1Len-k:f1Len] == frag2[0:k]:
			contig = frag1[0:f1Len-k] + frag2
			overlap = k
			checker = True
			return contig,overlap,checker
		elif frag2[f2Len-k:f2Len] == frag1[0:k]:
			contig = frag2[0:f2Len-k] + frag1
			overlap = k
			checker = False
			return contig,overlap,checker
		k -= 1
	#end while
	return contig,overlap,checker

#Main
infile = open("file.txt",'r')

#Get user input for parameters
fragMin = int(input("Input minimum fragment size: "))
fragMax = int(input("Input maximum fragment size: "))
fold = int(input("Input minimum coverage fold: "))

title = infile.readline()#Read first line of Fasta formatted filed
sequence = ""
for line in infile:#Read in remaining data in Fasta file
	line.replace("\n","")
	line.replace("\r","")
	sequence = sequence + line

coverage=[0] * len(sequence) #Holds coverage count of nucleotides
numFrags = 0
frags = [] #Create list of fragments

while(not coverageMet(coverage,fold)): #Continue until the coverage is met
	randLength = randint(fragMin,fragMax) 
	randStart = randint(0,len(sequence)-randLength)
	newFrag = (sequence[randStart:randStart+randLength])#Creates random fragment of the sequence 
	frags.append(newFrag) 
	numFrags +=1
	for i in range(randStart,randStart+randLength):
		coverage[i] += 1 #Increase cover 
#end while

newFragmentList = list() #List to check for substrings

#Ensure the new frag is not a substring of an old frag and an old frag is not a substring of the new frag
for fragCheck in frags:
	if all(fragCheck not in frag for frag in newFragmentList) and all(frag not in fragCheck for frag in newFragmentList):
		newFragmentList.append(fragCheck)

frags = newFragmentList#Set list with all substrings removed back to frags

fragDict = createDict(frags)
prefix = ()
suffix = ()
suffixList = list()
prefixList = list()

while len(prefixList)+1 < len(frags):
	maxValue = findBiggest(fragDict,frags)
	if maxValue == 0:#Makes sure there are sequences that overlap at one point
		break
	#finds biggest overlap in the list
	keys = findKeys(fragDict,maxValue)
	prefix = keys[0]
	suffix = keys[1]

	prefixList.append(prefix)
	suffixList.append(suffix)
	for k,v in fragDict.items():
		for k2,v2 in v.items():
			if k2 == (suffix):
				del fragDict[k][k2]
	del fragDict[prefix]
#end while

#Find what is in prefixlist and not in suffix list to start the sequence 
firstNode = [item for item in prefixList if item not in suffixList]
print("Fragments for sequence...")
print(frags)
sequenceList = [0] * len(firstNode)
for i in range(0,len(firstNode)): #Creates multiple fragments if no optimal solution exists 
	start = firstNode[i] 

	suffixPointer = prefixList[0]
	prefixPointer = start
	fullSequence = frags[start]

	while suffixPointer in prefixList: #Make sure there is a pointer to the next prefix if not end
		suffixPointer = suffixList[prefixList.index(prefixPointer)]
		suffix = frags[suffixPointer]
		fullSequence,overlap,checker = findOverlap(fullSequence,suffix)
		prefixPointer = suffixPointer
	sequenceList[i] = fullSequence
	#end while

if len(firstNode) == 1: #Make sure there is 1 starting position
	print(title + "Aligned sequence:")
	print(fullSequence)
else:
	print(title + "No optimal solution fragments of:") 
	for i in range(0,len(firstNode)):
		print("Fragment " + str(i +1) + " " + sequenceList[i])




