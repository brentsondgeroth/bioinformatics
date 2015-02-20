'''
Brent Gaither Gerardo Paleo

This program creates fragments of DNA from a Fasta formatted file. The fragments
are then overlapped to recreate the original fragment. 
Due Date 2/24/15
'''
from random import randint
from collections import defaultdict
import operator

def calculate_dict_max(dict_name):
    result = max(zip(dict_name.values(),dict_name.keys()))
    return result

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

def coverageMet(coverage,fold):
	i = 0
	met = True
	while(i<len(coverage) and met == True):
		if coverage[i] < fold:
			met = False
		i += 1
	return met

#Main

infile = open("file.txt",'r')

#Get user input for parameters
fragMin = int(raw_input("please input minimum fragment size "))
fragMax = int(raw_input("please input maximum fragment size "))
fold = int(raw_input("please input fold count "))

infile.readline()#Read first line of Fasta formatted filed
sequence = ""
for line in infile:
	line.replace("\n","")
	line.replace("\r","")
	sequence = sequence + line

coverage=[0] * (len(sequence) -1)

numFrag = 0
frag = list()#Create list of fragments
fragDict = defaultdict(lambda:defaultdict(int))

while(coverageMet(coverage,fold) is False):#Continue while until the coverage is high enough
	randLength = randint(fragMin,fragMax)
	randStart = randint(0,len(sequence)-randLength)

	newFrag = (sequence[randStart:randStart+randLength])
	if newFrag not in frag:
		frag.append(newFrag) 
		numFrag +=1
		for i in range(randStart,randStart+randLength-1):
			coverage[i] += 1
	else:
		print("Frag already contained")
		print(newFrag)
print(frag)

for i in range(0,numFrag-1):
	for j in range(i+1,numFrag):
		f1Len = len(frag[i])
		f2Len = len(frag[j])
		minLen = min(f1Len,f2Len)
		overlap = 0
		frag1 = frag[i]
		frag2 = frag[j]
		k = minLen -1 

		while k >= 1 and overlap == 0: #check for overlap of sequences and print out all overlapping sequences

			if frag1[f1Len-k:f1Len] == frag2[0:k]:
				contig = frag1[0:f1Len-k] + frag2
				overlap = k
				fragDict[frag1][frag2] = overlap
				print(frag1)
				print(frag2)
				print(contig)
				print(overlap)
			elif frag2[f2Len-k:f2Len] == frag1[0:k]:
				contig = frag2[0:f2Len-k] + frag1
				overlap = k
				fragDict[frag2][frag1] = overlap
				print(frag1)
				print(frag2)
				print(contig)
				print(overlap)
			k -= 1
print(fragDict)

prefix = list()
suffix = list()
#NEED TO FIND MAX IN DICT
biggest = calculate_dict_max(fragDict)

#key = findKeys(fragDict,max(fragDict))
print("BIGGEST")
print(biggest)

#print(max(fragDict))
#print(key)
