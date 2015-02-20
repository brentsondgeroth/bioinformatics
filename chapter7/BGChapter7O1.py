'''
Names:  Brent Gaither
Program that merges clusters using neighbor joining and outputs the final
product in Newick format with branch lengths.

Due Date:  02/19/2015
'''

from __future__ import division
import collections
from collections import defaultdict
from collections import Mapping
import math
import sys

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
readIn- Reads in a fasta formatted file into a list 
'''

def readIn(infile):
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

    smallest = sys.maxint
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

# main

# Read in data for nested hash structure
outfile = open("bgwmCh7skillsout.txt", 'w')
infile = open("file3.txt",'r')
sequenceList,clusterNames = readIn(infile)
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


