'''
Brent Gaither, Gerardo Paleo

This program creates fragments of DNA from a Fasta formatted file. The fragments
are then overlapped to recreate the original fragment. 
Due Date 2/24/15
'''
from random import randint

#Function to determine if coverage has been met
def coverageMet(coverage,fold):
	i = 0
	met = True #Assume coverage met
	while(i<len(coverage) and met):
		if coverage[i] < fold:
			met = False
		i += 1
	return met

#Prompt user for filename containing sequence and open it
filename = raw_input("Enter .txt filename that contains sequence: ")
infile = open(filename,'r')

#Get user input for parameters
fragMin = int(raw_input("Input minimum fragment size: "))
fragMax = int(raw_input("Input maximum fragment size: "))
fold = int(raw_input("Input minimum coverage fold: "))

infile.readline() #Discard first line of Fasta formatted file
sequence = ""
for line in infile:
	line.replace("\n","")
	line.replace("\r","")
	sequence = sequence + line

coverage=[0] * len(sequence) #Holds coverage count of nucleotides

numFrags = 0
frags = [] #Holds list of fragments

#Generate a set of fragments for the input sequence. Stop once coverage is met. 
while(not coverageMet(coverage,fold)): 
	randLength = randint(fragMin,fragMax)
	randStart = randint(0,len(sequence)-randLength)
	newFrag = (sequence[randStart:randStart+randLength])
	#Accept new fragment if it is not a substring of a fragment already on the list or viceversa
	if all(newFrag not in frag for frag in frags) and all(frag not in newFrag for frag in frags):
		frags.append(newFrag) 
		numFrags +=1
		#Update coverage
		for i in range(randStart,randStart+randLength):
			coverage[i] += 1	
print "\nAccepted Fragments:"
print(frags)

#Determine overalp for each pair of fragments
for i in range(0,numFrags-1):
	for j in range(i+1,numFrags):
		f1Len = len(frags[i])
		f2Len = len(frags[j])
		minLen = min(f1Len,f2Len)
		overlap = 0
		frag1 = frags[i]
		frag2 = frags[j]
		k = minLen -1 
		while k >= 1 and overlap == 0: 
			#Compare suffix of frag1 to prefix of frag2
			if frag1[f1Len-k:f1Len] == frag2[0:k]:
				#Create contig
				contig = frag1[0:f1Len-k] + frag2
				overlap = k
				print "\nFragment 1: " + (frag1)
				print "Fragment 2: " + (frag2)
				print "Contig:     " + (contig)
				print "Overlap:    " + str(overlap)
			#Compare suffix of frag2 to prefix of frag1 
			elif frag2[f2Len-k:f2Len] == frag1[0:k]:
				#Create contig
				contig = frag2[0:f2Len-k] + frag1
				overlap = k
				print "\nFragment 1: " + (frag1)
				print "Fragment 2: " + (frag2)
				print "Contig:     " + (contig)
				print "Overlap:    " + str(overlap) 
			k -= 1