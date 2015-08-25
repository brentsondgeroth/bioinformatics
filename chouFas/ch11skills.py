'''
Ch 11 Skills 1 - 2
by Leiyang G.  and Brent G.
This program take in a sequence and output an aligned string of alpha helix
'''
# Step 1: read in sequence from a file
aminoSeq = "" 
infile = open(raw_input("Enter the file name of the sequence (no suffix and must be in FASTA format):") + ".txt","r") 
infile.readline()  # bypass > header line
for line in infile:
	line = line.replace('\n', '')
	aminoSeq = aminoSeq + line
aminoSeq = aminoSeq.upper()

# Step 2: find Alpha Helices
# find region of six (step 1a)
# initializing all values and hashes
lenSeq = len(aminoSeq)
window = 6
pScore = 103
minWindow = 4
paHash = {"A":142, "R":98, "N":67, "D":101, "C":70, "E":151, "Q":111, "G":57, "H":100, "I":108, "L":121, "K":114, "M": 145, "F":113, "P":57, "S":77, "T":83, "W":108, "Y":69, "V":106}
pbHash = {"A":83, "R":93, "N":89, "D":54, "C":119, "E":37, "Q":110, "G":75, "H":87, "I":160, "L":130, "K":74, "M": 105, "F":138, "P":55, "S":75, "T":119, "W":137, "Y":147, "V":170}
helixStart = []
helixEnd = []
newStart = 0

for i in range(lenSeq - window + 1):
	# check for helix overlapping (to ensure no 2 regions will overlap)
	if i > newStart:
		ctr = paSum = pbSum = 0
		# calculating sum of P(a) and sum of P(b) values in a region
		for j in range(window):
			paSum = paSum + paHash[aminoSeq[i+j]]
			pbSum = pbSum + pbHash[aminoSeq[i+j]]
			if paHash[aminoSeq[i+j]] > pScore:
				ctr += 1
		
		# check if a possible region found
		if ctr >= minWindow:
			extend = i - 1
			done = False
			
			# check if the region found can be an actual alpha helix region
			while extend >= 0 and done == False:
				if extend >= 3 and paHash[aminoSeq[extend]] < 100 and paHash[aminoSeq[extend-1]] < 100 and paHash[aminoSeq[extend-2]] < 100 and paHash[aminoSeq[extend-3]] < 100:
					done = True
				else:
					paSum = paSum + paHash[aminoSeq[extend]]
					pbSum = pbSum + pbHash[aminoSeq[extend]]
					extend -= 1
			left = extend + 1
			extend = i + window
			done = False
			
			while extend < lenSeq and done == False:
				if extend <= lenSeq - 3 and paHash[aminoSeq[extend]] < 100 and paHash[aminoSeq[extend+1]] < 100 and paHash[aminoSeq[extend+2]] < 100 and paHash[aminoSeq[extend+3]] < 100:
					done = True
				else:
					paSum = paSum + paHash[aminoSeq[extend]]
					pbSum = pbSum + pbHash[aminoSeq[extend]]
					extend += 1
			right = extend - 1
			
			# last step of checking: if the sum is greater than threshold value
			lenRegion = right - left
			if paSum/lenRegion > pScore and paSum > pbSum:
				helixStart.append(left)
				helixEnd.append(right)
				newStart = right + 1

dirString = ""

# Part 2 Skills 2 building helix string and align it with amni sequence in output
for i in range(lenSeq):
	dirString = dirString + "-"

	for j in range(len(helixStart)):
		if i >= helixStart[j] and i <= helixEnd[j]:
			dirString = dirString[:len(dirString) - 1]
			dirString = dirString + "H"
			break

print "Pred:   " + dirString
print "AA:     " + aminoSeq

			
		
			
