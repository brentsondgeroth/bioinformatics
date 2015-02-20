'''
Name: Brent Gaither Kevin Portland
Description: This program uses the needleman wunch algorithm to align amino acid
sequences globally using a hashmap for amino acid similarities. 
bgkpCh5S1
Due Date: 1/29/15
'''

import collections
from collections import defaultdict

'''
printMatrix- Uses the lists that make up the matrix and formats
to print in a table.
'''
def printMatrix(matrix,outfile):
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print ('\n'.join(table))
    outfile.write('\n'.join(table))
'''
buildMatrix- Builds a string that provides direction for tracing
a path back through the matrix to build the best match.
'''
def buildMatrix(matrix,gap,N,M,s1,s2,s):

    #Step 1 build matrix
    matrix [0][0] = 0
    for i in range(1,N+1): 
        matrix[i][0] = matrix[i-1][0]+float(gap) #Sets col to penalize for initial or terminal gaps
        for j in range(1,M+1):
            score1 = matrix[i - 1][j - 1] + s[s1[i - 1]][s2[j - 1]]  #Use dict to find amino acid similarity 
            score2 = matrix[i][j - 1] + float(gap)
            score3 = matrix[i - 1][j] + float(gap)
            matrix[i][j] = max(score1, score2, score3) 

    return matrix
'''
buildDirectional- Creates the directional string for following the matrix to 
the best solution
'''
def buildDirectional(matrix, N, M, gapScore):
    dString = ''
    currentRow = N
    currentCol = M
    while(currentRow != 0 or currentCol != 0):

        if(currentRow == 0):
            dString = dString + ('H')
            currentCol = currentCol - 1
        elif(currentCol == 0):
            dString = dString + ('V')
            currentRow = currentRow - 1
        elif(matrix[currentRow][currentCol - 1] + float(gapScore) == matrix[currentRow][currentCol]):
            dString = dString + ('H')
            currentCol = currentCol - 1
        elif(matrix[currentRow - 1][currentCol] + float(gapScore) == matrix[currentRow][currentCol]):
            dString = dString + ('V')
            currentRow = currentRow - 1
        else:
            dString = dString + ('D')
            currentRow = currentRow - 1
            currentCol = currentCol - 1
    return dString
'''
readIn- Reads in the files
'''
def readIn(infile):
    sequence = ""
    infile.readline()  # bypass > header line
    for line in infile:
        line = line.replace('\n', '')
        sequence = sequence + line
    sequence = sequence.upper()
    return sequence
'''
Main
'''
# Opens files
infile1 = open('yadk.txt', 'r')
infile2 = open('yadk2.txt', 'r')
outfile = open('ch5sk1out.txt', 'w')

#Prompt use for score preferences
print("This program aligns amino acid sequences.")
pickMatrix = input ("Which substitution matrix would you like to use: BLOSUM62(1) PAM250(2) PAM1(3)\n" +
 "Hydrophobicity(4) User created substitution matrix(5) ")
#open substitution matrix
if pickMatrix == '1':
    infile3 = open('BLOSUM62.txt','r')
elif pickMatrix == '2':
    infile3 = open('PAM250.txt','r')
elif pickMatrix == '3':
    infile3 = open('PAM1.txt','r')
elif pickMatrix == '4':
    infile3 = open('Hydrophobicity.txt','r')
else:
    infile3 = open('ch5O1out.txt','r')

gap = input("Please enter gap score: ")

# Reads the files into strings.
s1 = readIn(infile1)
s2 = readIn(infile2)

temp = ""
if (len(s1) > len(s2)): #ensure long sequence is s2
    temp = s1
    s1 = s2
    s2 = temp

#Initialize matrix
N = len(s1)
M = len(s2)
connectors = ''
matrix = ''
matrix = list(matrix)
matrix = [[0 for col in range(M+1)] for row in range(N+1)]

s = defaultdict(lambda : defaultdict(dict)) #Create anynomous function

aminoAcids = ''
aminoAcidArray = []
infile3.readline() # bypass > header line
aminoAcids = infile3.readline() #get amino acid order
aminoAcidArray = aminoAcids.strip().split(',') 

i = 0
for line in infile3: #Read in the substitution matrix place in dict
    j = 0
    lineList = list(line.strip().split(","))
    for value in lineList:
        s[aminoAcidArray[i]][aminoAcidArray[j]] = float(value)
        j += 1
    i += 1

infile3.close()

#Build matrix with input sequences
matrix = buildMatrix(matrix,gap,N,M,s1,s2,s)

#step 2 create directional string
dString = ''
dString = buildDirectional(matrix, N, M, gap)
#step 3 Build alignment using directional strings
seq1Pos = N-1
seq2Pos = M-1
dirPos = 0
alSeq1 = ''
alSeq2 = ''
alignmentScore = 0
matrixOut = ''

while(dirPos < len(dString)):
        
    if(dString[dirPos] == "D"): #Align sequence for match
        alSeq1 = s1[seq1Pos] + alSeq1
        alSeq2 = s2[seq2Pos] + alSeq2
        seq1Pos = seq1Pos - 1
        seq2Pos = seq2Pos - 1
    elif(dString[dirPos] == "V"):#Place a space for a gap
        alSeq1 = s1[seq1Pos] + alSeq1
        alSeq2 = '-' + alSeq2
        seq1Pos = seq1Pos - 1
    else:
        alSeq1 = '-' + alSeq1 #Place a space for a gap
        alSeq2 = s2[seq2Pos] + alSeq2
        seq2Pos = seq2Pos -1

    dirPos = dirPos + 1

alSeq1List = list(alSeq1)
alSeq2List = list(alSeq2)

for char1, char2 in zip(alSeq1List, alSeq2List):
    if char1 == char2:
        connectors = connectors + "|" #Create list for relationship between sequences
        alignmentScore = alignmentScore + 1
    else:
        if char1 == "-" or char2 == "-":
            connectors = connectors + " "
        else:
            connectors = connectors + "."

line = 0
count = 1
charPerLine = 50

# Prints the alignment with a limit of 50 characters per line.
for line in range(line,len(alSeq1),charPerLine):
    print(str(count) + " " + alSeq1[line:line + charPerLine])
    outfile.write(str(count) + " " + alSeq1[line:line + charPerLine]+ '\n')
    print(len(str(count)) * " " + " " + connectors[line:line + charPerLine])
    outfile.write(len(str(count)) * " " + " " + connectors[line:line + charPerLine]+ '\n')
    print(str(count) + " " + alSeq2[line:line+charPerLine])
    outfile.write(str(count) + " " + alSeq2[line:line+charPerLine])
    count = count + charPerLine

print("Match percentage is: " + str(alignmentScore/float(len(alSeq1))*100))
outfile.write("\nMatch percentage is: " + str(alignmentScore/float(len(alSeq1))*100) + '\n')
print("The alignment score is: " + str(alignmentScore))
outfile.write("The alignment score is: " + str(alignmentScore)+ '\n' )
matrixOut = input("Would you like to print the matrix? (1 yes, 2 no) ")
if (matrixOut == "1"):
	printMatrix(matrix,outfile)
dStringOut = input("Would you like to print the path string? (1 yes, 2 no) ")
if (dStringOut == "1"):
    print (dString)
    outfile.write('\n' + dString)

