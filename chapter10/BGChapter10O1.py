'''
Brent Gaither
This program uses a Hidden Markov model to calculate the probability of a 
eukaryotic sequence. 
HMM 

Due Date 3/10/15
'''

from collections import defaultdict
pathList = list()
p = list()
def readIn (infile):
	infile.readline()#Read first line of Fasta formatted filed
	sequence = ""
	for line in infile:#Read in remaining data in Fasta file
		line.replace("\n","")
		line.replace("\r","")
		sequence = sequence + line
	sequence = sequence.upper()
	return sequence

def markov (states,currentState, length, ctr,end, path,sequence,emmision,transition):
	#if currentState	== end or ctr == length:
	if ctr >= length:
		if len(path) > len(sequence):
			return
		'''
		prob = dict()
		for state in path:
			print(sequence[ctr])
			prob[state] = float(transition[state][currentState]) * float(emmision[state][sequence[ctr-1]])
		print((prob))
		'''
		'''
		for state in path:
			p[i] = prob[state]
		'''
		'''
			print(path)
			print("trans")
			print(state)
			print(currentState)
			print(transition[state][currentState])
			print("ems")
			print(emmision[state][sequence[ctr-1]])
		'''
			#print((float(transition[state][currentState]) * float(emmision[state][sequence[ctr-1]])))
		return pathList.append(path)
	else:
		for state in states[currentState]:	
			path = path + state
			markov(states,state, len(sequence), len(path),'X', path,sequence,emmision,transition)

#main
infile = open("file1.txt",'r')

sequence = readIn(infile)

transition = defaultdict(lambda:defaultdict(float))
emmision = defaultdict(lambda:defaultdict(float))

states = defaultdict(lambda:defaultdict())
'''
B = begin, A = alpha Exon, E = internal exon, S,P splice donor sites, 
C,R splice acceptor site, I = internal exon FGH = first ATG, XYZ = terminal nucleotides of terminal exon
'''
#Begin State 
states['B'] = ['A']
#Alpha Exon state
states['A'] = ['A', 'F']
#First ATG state
#state = {'F':'G','G':'H','H':'E'}
states['F'] = ['G']
states['G'] = ['H']
states['H'] = ['E']

#Splice state
states['S'] = ['P']
states['P'] = ['I']
states['C'] = ['R']
states['R'] = ['E']

#Exon states
states['E'] = ['E','S','X']
states['I'] = ['I','C']
'''
states['X'] = ['Y']
states['Y'] = ['Z']
states['Z'] = ['']
'''
ctr = 0

transition['E']['E'] = .8
transition['E']['S'] = .1
transition['E']['X'] = .1
transition['I']['I'] = .9
transition['S','P'] = 1
transition['I']['S'] = .1
transition['P']['E'] = 1
transition['A']['A'] = .9
transition['A']['F'] = .1
transition['F']['G'] = 1
transition['G']['H'] = 1
transition['H']['E'] = 1



emmision['C']['A'] = .9998
emmision['C']['T'] = .000067
emmision['C']['G'] = .000067
emmision['C']['C'] = .000067
emmision['R']['A'] = .0005
emmision['R']['T'] = .0001
emmision['R']['G'] = .9993
emmision['R']['C'] = .0001

emmision['S']['A'] = .0005
emmision['S']['T'] = .0001
emmision['S']['G'] = .9993
emmision['S']['C'] = .0001
emmision['P']['A'] = .0001
emmision['P']['T'] = .0069
emmision['P']['G'] = .0001
emmision['P']['C'] = .9929

emmision['F']['A'] = 1
emmision['G']['T'] = 1
emmision['H']['G'] = 1

emmision['E']['A'] = .2
emmision['E']['T'] = .3
emmision['E']['G'] = .3
emmision['E']['C'] = .2

emmision['A']['A'] = .2
emmision['A']['T'] = .3
emmision['A']['G'] = .3
emmision['A']['C'] = .2

markov(states,'A',len(sequence),ctr,'X','',sequence,emmision,transition)
'''
for i in range (0,len(p))
	p[i] = p[i] * float(transition[state][currentState]) * float(emmision[state][sequence[ctr-1]])
'''

j = 0
print(len(sequence))
total = [float(1)] *(len(pathList))
print(total)
print(pathList)
for path in pathList:
	previousState = 'A'
	i = 0
	for state in path:
		'''
		print(state)
		print(previousState)
		print(sequence[i])
		print(float(emmision[state][sequence[i]]) * float(transition[previousState][state]))
		'''
		p.append(float(transition[previousState][state]) * float(emmision[state][sequence[i]]))
		previousState = state
		i+=1
	for num in p:
		total[j] *= float(num)
	j+=1
	del p[:]
	print("NEW")


'''
for j in range(0,len(sequence)):
	for num in p:
		print( total[j] * num)
		total[j] = float(total[j]) * float(num)
'''
print(total)
print(max(total))
print("Most likely start eukaryotic sequence is...")
print(pathList[(total.index(max(total)))])
#print(sequence[pathList[(total.index(max(total))).index("F")]])

#print(max(p))


