'''
Names:  Brent Gaither
Program that merges clusters using neighbor joining and outputs the final
product in Newick format.

Due Date:  02/19/2015
'''

from __future__ import division
from collections import defaultdict
from collections import Mapping
import math
import sys

'''
compareSequences- reads in the sequences and compares to find matches, length of sequences,
mutation type and the GC content. 
'''

def compareSequences(alSeq1,alSeq2,choice):
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

    if choice == '1':
        return findJukesCantor(difCtr,lenCtr)
    elif choice == '2':
        return findKimura(S,V,lenCtr)
    else:
        return findTamura(S,V,lenCtr,gcContent1,gcContent2)
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
countSequences- Counts number fasta formatted sequences in the infile
'''

def countSequences(infile):

    in_FH = open('file3.txt','r')
    counter=0
    for line in in_FH:
        if line[:1] is '>':
            counter += 1
    return counter

'''
readIn- Reads in fasta formatted file into a list 
'''

def readIn(infile):
    sequenceList = [0] * countSequences(infile)
    sequence = ''
    i = 0
    first = True
    for line in infile:
        if line[0]== '>' and first != True:
            print(line)
            print(sequence)
            sequenceList[i] = sequence
            sequence = ''
            i += 1
        else:
            line = line.replace('\n', '')
            line = line.replace('\r', '')
            sequence = sequence + line
            first = False
        sequence = sequence.upper()
    return sequenceList
# Builds a reverse copy of a nested dictionary
def keypaths(nested):
    for key, value in nested.iteritems():
        if isinstance(value, Mapping):
            for subkey, subvalue in keypaths(value):
                yield[key] + subkey, subvalue
        else:
            yield[key], value

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

def findSmallest(clusterDist):
    smallest = sys.maxint
    for k,v in clusterDist.items():
        for k2,v2 in v.items():
            if k != k2:
                if smallest > int(v2) and int(v2) != 0:
                    smallest = int(v2)
    return smallest
       
'''
singleLinkage- Calculate distances between new cluster and all
other clusters using single linkage
'''

def singleLinkage(clusterDist, newCluster, originalDist, clusterNames, reverseCluster):
    for cluster in clusterNames:
        smallestD = max(reverseCluster)
        for c1 in cluster:
            for c2 in newCluster:
                if originalDist[c1][c2] < smallestD:
                    smallestD = originalDist[c1][c2]
        clusterDist[newCluster][cluster] = smallestD
        clusterDist[cluster][newCluster] = smallestD
    return

'''
calculateTransition- Calculates the transition matrix for a neighbor joining phylogeny.
'''

def calculateTransition(clusterDist,transition):
    r = [0] * len(clusterNames)
    dix = 0
    for i in range (0, int(numSeq)):
        for j in range (0, int(numSeq)):
            dix = dix + int(clusterDist[clusterNames[i]][clusterNames[j]])
            j += 1
        r[i] = dix / (len(clusterNames) - 2)
        i +=1
        dix  = 0
    for i in range (0, int(numSeq)):
        for j in range (0, int(numSeq)):
            transition[clusterNames[i]][clusterNames[j]] = int(clusterDist[clusterNames[i]][clusterNames[j]])  - r[i] - r[j]

def neighborJoining(clusterDist, newCluster, originalDist, clusterNames, reverseCluster):
    for cluster in clusterNames:
        dk = (originalDist[clusterNames[shortestI]][clusterNames[shortestD]] + originalDist[clusterNames[shortestJ]][clusterNames[shortestD]] - originalDist[clusterNames[shortestI]][clusterNames[shortestJ]]) /2
# main

# Read in data for nested hash structure
infile1 = open("file2.txt", 'r')
#infile2 = open('file3.txt', 'r')
outfile = open("bgwmCh7skillsout.txt", 'w')

# sequenceList = readIn(infile2)

numSeq = infile1.readline()
clusterNames = [0] * int(numSeq)
distances = [0] * int(numSeq)
idx = 0
for line in infile1:
    clusterNames[idx] = line[0]
    distances[idx] = line[1:len(line) - 1:1].split()
    idx += 1

# Build nested hash structures
originalDist = defaultdict(lambda:defaultdict(int))
clusterDist = defaultdict(lambda:defaultdict(int))
transition = defaultdict(lambda:defaultdict(int))

for i in range (0, int(numSeq)):
    for j in range (0, int(numSeq)):
        originalDist[clusterNames[i]][clusterNames[j]] = distances[i][j]
        clusterDist[clusterNames[i]][clusterNames[j]] = distances[i][j]
# Step 1:  Cluster

# Set flag equal to the total number of clusters
numClusters = len(clusterNames)

'''
Create a reverse copy of clusterDist to be passed with
the shortest distance into findKeys to return the parent
and child keys of the shortest distance
'''
reverseCluster = {}
for keypath, value in keypaths(clusterDist):
    reverseCluster.setdefault(value, []).append(keypath)

# Initialize a nested dictionary to build the Newick format output
newick = defaultdict(lambda:defaultdict(float))

calculateTransition(clusterDist,transition)

while(numClusters > 2):
    shortestD = findSmallest(clusterDist)
    # Finds keys of smallest value
    shortestK = findKeys(clusterDist, str(shortestD))
    shortestK = sorted(shortestK)
    shortestI = shortestK[0]
    shortestJ = shortestK[1]
    print(shortestJ)

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
    myValue = '(' + x + ',' + y + ')'

    newick[myKey][myValue] = 0
    
    # merge clusters I and J
    newCluster = shortestI + shortestJ
    clusterNames.remove(shortestI)
    clusterNames.remove(shortestJ)

    # First, remove the child key and values
    for k,v in clusterDist.items():
        for k2,v2 in v.items():
            if k2 == shortestI:
                del clusterDist[k][k2]
            if k2 == shortestJ:
                del clusterDist[k][k2]

    # Next, merge the 2 clusters
    #singleLinkage(clusterDist, newCluster, originalDist, clusterNames, reverseCluster)
    neighborJoining(clusterDist, newCluster, originalDist, clusterNames, reverseCluster)
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

    print "merging clusters " + shortestI + " and " + shortestJ
    outfile.write("merging clusters " + shortestI + " and " + shortestJ + '\n')

    numClusters -= 1

# Put child keys from Newick in list, then join list elements as a string
nFormat = []

for k,v in newick.items():
    for k2,v2 in v.items():
        nFormat.append(k2)

nFormat = ','.join(nFormat)


print sorted(clusterNames)

# Surround the newick string with a final set of parentheses
print "(" + nFormat + ")"
outfile.write(','.join(clusterNames))
outfile.write("\n(" + nFormat + ")")
outfile.close()


