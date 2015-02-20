'''
Name: Matthew Kachlik
Due: 2/10/15
Description: Implements Jukes-Cantor model 
'''
import math

'''
	Jukes
	measures evolutionary distance between two nucleotide sequences
'''
def Juke(seq1,seq2):

	difCtr = float(0)
	lenCtr = float(0)

	#count differences between sequences ignoring gaps
	for i in range(0, len(seq2)):
		if(seq1[i] != '-' and seq2[i] != '-'):
			lenCtr += 1
			if(seq1[i] != seq2[i]):
				difCtr += 1
	print(difCtr)
	print(lenCtr)
	#Calculate and output results
	difCtr = difCtr/lenCtr
	print(difCtr)
	jukes = (-3/4) * math.log(1-(4/3*difCtr))
	print (jukes)
	return jukes


'''
	main
'''

# Opens files
infile2 = open('giraffe.txt', 'r')
infile1 = open('hippo.txt', 'r')

# Reads the files into strings.
s1 = ""
infile1.readline()  # bypass > header line
for line in infile1:
    line = line.replace('\n', '')
    s1 = s1 + line
s1 = s1.upper()
infile1.close()

s2 = ""
infile2.readline()  # bypass > header line
for line in infile2:
    line = line.replace('\n', '')
    s2 = s2 + line
s2 = s2.upper()
infile2.close()
distance = Juke(s1,s2)
