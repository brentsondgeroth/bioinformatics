'''
Brent Gaither Lei Guo
'''


'''
readIn- Method used to read in a FASTA formatted file
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

#main
infile = open("1KJFsequence.txt",'r')
chouFasFile = open("chouFas.txt",'r')

acids = list()
paScore = list()
pbScore = list()

chouFas = chouFasFile.read().split()

for i in range (0,len(chouFas),3):
	acids.append(chouFas[i]) 
	paScore.append(chouFas[i+1]) 
	pbScore.append(chouFas[i+2]) 
sequence = readIn(infile)

lenSeq = len(sequence)
window = 6
pScore = 103
minWindow = 4

paHash={}
pbHash={}

for i in range(len(acids)):
    paHash[acids[i]] = paScore[i]
    pbHash[acids[i]] = pbScore[i]

for i in range(0,lenSeq-window):
	ctr = 0
	paSum = 0
	pbSum = 0
	for j in range(0,window-1):
		paSum = paSum + int(paHash[sequence[i+j]])
		pbSum = pbSum + int(pbHash[sequence[i+j]])
		if int(paHash[sequence[i+j]]) > pScore:
			ctr +=1
	if ctr>= minWindow:
		print("Possible alpha helix region found at " + str(i+1))
		extend = i-1
		done = False

		while extend >= 0 and not done:
			if extend >=3 and int(paHash[sequence[extend]])<100 and int(paHash[sequence[extend-1]])<100 and int(paHash[sequence[extend-2]])<100 and int(paHash[sequence[extend-3]])<100:
				done = True
			else:
				paSum = paSum + int(paHash[sequence[extend]])
				pbSum = pbSum + int(pbHash[sequence[extend]])
				extend -=1
		left = extend +1
		extend = i + window
		done = False
		while extend<lenSeq and not done:
			if extend <= lenSeq-3 and int(paHash[sequence[extend]])<100 and int(paHash[sequence[extend-1]])<100 and int(paHash[sequence[extend-2]])<100 and int(paHash[sequence[extend-3]])<100:
				done = True
			else:
				paSum = paSum + int(paHash[sequence[extend]])
				pbSum = pbSum + int(pbHash[sequence[extend]])
				extend +=1
		right = extend -1
		lenRegion = right - left
		if (paSum/lenRegion) > pScore and paSum > pbSum:
			print("Alpha region: " + (leftStart+1) + " to " + (rightStart+1))







