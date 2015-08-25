'''
Names:  Brent Gaither and Bill Moers
Program that merges clusters by distance and outputs the final
product in Newick format.

bgwmCh7S2

Due Date:  02/17/2015
'''

from __future__ import division
from collections import defaultdict
import math
import sys

'''
Used the smallest value and the reverse copy of the nested dictionary
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
Finds and returns the shortest distance between sequences
listed in the matrix
'''
def findSmallest(clusterDist):
    smallest = sys.maxint
    for k,v in clusterDist.items():
        for k2,v2 in v.items():
            if k != k2:
                if smallest > int(v2) and int(v2) != 0:
                    print("HERE")
                    smallest = int(v2)

    return smallest
       
'''
Calculate distances between new cluster and all
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
    return

# main

# Read in data for nested hash structure
infile1 = open("file2.txt", 'r')
outfile = open("bgwmCh7skillsout.txt", 'w')
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
for i in range (0, int(numSeq)):
    for j in range (0, int(numSeq)):
        originalDist[clusterNames[i]][clusterNames[j]] = distances[i][j]
        clusterDist[clusterNames[i]][clusterNames[j]] = distances[i][j]

# Step 1:  Cluster

# Set flag equal to the total number of clusters
numClusters = len(clusterNames)

# Initialize a nested dictionary to build the Newick format output
newick = defaultdict(lambda:defaultdict(float))

while(numClusters > 2):
    shortestD = findSmallest(clusterDist)
    
    # Finds keys of smallest value
    shortestK = findKeys(clusterDist, str(shortestD))
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
    singleLinkage(clusterDist, newCluster, originalDist, clusterNames)
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


print (sorted(clusterNames))

# Surround the newick string with a final set of parentheses
print ("(" + nFormat + ")")
outfile.write(','.join(clusterNames))
outfile.write("\n(" + nFormat + ")")
outfile.close()

