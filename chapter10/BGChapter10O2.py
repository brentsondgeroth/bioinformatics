'''
Brent Gaither
This program uses a Hidden Markov model to calculate the probability of a 
eukaryotic sequence. 

Due Date 3/10/15
'''

from collections import defaultdict
import math

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

'''
markov- recursive method to create all the possible paths for 
the sequence
'''

def markov (states,currentState, length, ctr,end,path):

	if ctr == length and currentState == end and path not in pathList:
		pathList.append(path)
	elif currentState == end:
		return 
	elif ctr<length:
		for state in states[currentState]:#Go through each possible state to create each unique path		
			markov(states,state, length, ctr+1,end,path+state)
	else:
		return
#main
infile = open("file1.txt",'r')

sequence = readIn(infile)

#Create nested hash tables
transition = defaultdict(lambda:defaultdict(float))
emission = defaultdict(lambda:defaultdict(float))
states = defaultdict(lambda:defaultdict())

'''
B = begin, W = start of sequence, A = alpha Exon, E = internal exon, S,P splice donor sites,
C,R splice acceptor site, I = internal intron FGH = first ATG, XYZ = terminal nucleotides of terminal exon
'''
states = {'B':['W'],'W':['W','F'],'A':['A','S'],'F':['G'],'G':['H'],'H':['A'],'E':['E','X','S'],'S':['P'],'P':['I'],'I':['I','C'],'C':['R'],'R':['E'],'X':['X']}

#Transition values
transition['W']['W'] = .9
transition['W']['F'] = .1
transition['F']['G'] = 1
transition['G']['H'] = 1
transition['H']['A'] = 1

transition['A']['A'] = .9
transition['A']['S'] = .1

transition['E']['E'] = .8
transition['E']['S'] = .1
transition['E']['X'] = .1
transition['S']['P'] = 1
transition['P']['I'] = 1

transition['I']['I'] = .9
transition['I']['C'] = .1
transition['C']['R'] = 1
transition['R']['E'] = 1

#emission values 
emission['C']['A'] = .9998
emission['C']['T'] = .000067
emission['C']['G'] = .000067
emission['C']['C'] = .000067
emission['R']['A'] = .0005
emission['R']['T'] = .0001
emission['R']['G'] = .9993
emission['R']['C'] = .0001

emission['S']['A'] = .0005
emission['S']['T'] = .0001
emission['S']['G'] = .9993
emission['S']['C'] = .0001
emission['P']['A'] = .0001
emission['P']['T'] = .0069
emission['P']['G'] = .0001
emission['P']['C'] = .9929

emission['F']['A'] = 1
emission['G']['T'] = 1
emission['H']['G'] = 1

emission['E']['A'] = .2
emission['E']['T'] = .3
emission['E']['G'] = .3
emission['E']['C'] = .2

emission['I']['A'] = .27
emission['I']['T'] = .3
emission['I']['G'] = .23
emission['I']['C'] = .2


emission['A']['A'] = .2
emission['A']['T'] = .3
emission['A']['G'] = .3
emission['A']['C'] = .2

emission['W']['A'] = .2
emission['W']['T'] = .3
emission['W']['G'] = .3
emission['W']['C'] = .2

emission['X']['A'] = 1
emission['X']['T'] = 1
emission['X']['G'] = 1
emission['X']['C'] = 1


ctr = 1
pathList = list()
total = list()
length = len(sequence)
#pathList = markov(states,'B',length,ctr,'X','')
markov(states,'B',length,ctr,'X','')

p = [float(1)] *(len(pathList))#Create list to hold probabilities 
j = 0
for path in pathList:
	previousState = 'W'#Set first previous state to sequence start
	i = 0
	for state in path:#multiply emission and transition values for total probability of a path
		total.append(float(transition[previousState][state]) * float(emission[state][sequence[i]]))
		previousState = state
		i+=1
	print(total)
	for num in total:
		p[j] *= float(num)
	j+=1
	del total[:]

k = 0
for num in p:#Calculate log of all probabilities
	if num!=0:
		p[k] = math.log(p[k])
	else:
		p[k] = 0
	k+=1

for i in range(0,len(pathList)):#Prints all possible paths
	print(pathList[i])

print("W = start of sequence, A = alpha Exon, E = internal exon, S,P splice donor sites," +
	"C,R splice acceptor site, I = internal intron FGH = first ATG, XYZ = terminal nucleotides of terminal exon")
print("Top sequence is the path list for the most likely eukaryotic start sequence (bottom)")

line = 0
nucleotideCount = 1
charPerLine = 50
startSite = max([x for x in p if x !=0])

# Prints the alignment with a limit of 50 characters per line.
for line in range(line,len(sequence),charPerLine):
    print(str(nucleotideCount) + " " + pathList[(p.index(startSite))][line:line + charPerLine])
    print(str(nucleotideCount) + " " + sequence[line:line+charPerLine])
    nucleotideCount = nucleotideCount + charPerLine
