'''
Name: Brent Gaither Bill Moers
Description: This program uses 
Due Date: 2/17/15
'''
import collections
from collections import defaultdict
'''
def singleLinkage():
	for cluster in clusterNames:
		smallestD = max int
		for c1 in cluster:
			for c2 in newClusterName:
				if originalDist[c1][c2]<smallestD:
		clusterDist[newClusterName][cluster] = smallestD
		clusterDist[cluster][newClusterName] = smallestD
'''

infile1 = open('samplephylip.txt', 'r')

numSeq = infile1.readline()
clusterNames = [0] * int(numSeq)
distances = [0] * int(numSeq)

i = 0
for line in range (0,int(numSeq)):
	first = infile1.readline()
	first = list(first)
	clusterNames[i] = first[0]
	distances[i] = first[1:len(first)]
	i +=1
'''
originalDist = defaultdict[int]
clusterDist = defaultdict[int]
'''
for i in range(0, int(numSeq)-1):
	for j in range (0,int(numSeq)-1):
		i +=1
		#originalDist[clusterNames[i]][clusterNames[j]] = distances[i][j]
	#clusterDist[clusterNames[i]][clusterNames[j]] = distances[i][j]

numClusters = #EQUAL TO?

while(numClusters>2):
	shortestD = min(clusterDist)
	shortestI = last(clusterDist)
	shortestJ = first(clusterDist)

	newClusterName = shortestI + shortestJ
	clusterNames.remove(shortestI)
	clusterNames.remove(shortestJ)
	del clusterDist[shortestI]
	del clusterDist[shortestJ]

	singleLinkage(clusterDist,newClusterName,originalDist,clusterNames)
	newClusterName.append(clusterNames) 

	print("merging clusters" + shortestI + shortestJ)

	numClusters -= 1

