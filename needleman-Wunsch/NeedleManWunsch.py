'''
Name: Brent Gaither 
Description: This program uses the needleman wunch algorithm to align 
sequences globally or semi globally or locally 
bgwmCh3O1
Due Date: 1/27/15
'''

def readInFile(infile):
    sequence = ""
    infile.readline()  # bypass > header line
    for line in infile:
        line = line.replace('\n', '')
        sequence = sequence + line
        sequence = sequence.upper()
    infile.close()
    return sequence

def checkSequenceLength(sequenceOne,sequenceTwo):

    if len(sequenceOne) > len(sequenceTwo):
        return sequenceOne,sequenceTwo
    else:
        transferSequence = sequenceOne
        sequenceOne = sequenceTwo
        sequenceTwo = transferSequence
        return sequenceOne, sequenceTwo

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
def buildMatrix(matrix, alignment, gap, misMatch, match, sequenceOne, sequenceTwo):

    #Step 1 build matrix
    N = len(sequenceOne)
    M = len(sequenceTwo)
    matrix [0][0] = 0
    for i in range(1,N+1): 
        if alignment == 'G':
            matrix[i][0] = matrix[i-1][0]+int(gap) #Sets col to penalize for initial or terminal gaps
        else:
            matrix [i][0] = matrix[i-1][0] #No penalty for local or semi local initial or terminal gaps
        for j in range(1,M+1):
            if alignment == 'G':
                matrix [0][j] = matrix[0][j - 1] + int(gap) 
            else:
                matrix [0][j] = matrix[0][j-1] 
            if (sequenceOne[i-1] == sequenceTwo[j - 1]):
                score1 = matrix[i - 1][j - 1] + int(match)               
            else:
                score1 = matrix[i - 1][j - 1] + int(misMatch)
            if alignment == 'G':
                score2 = matrix[i][j - 1] + int(gap)
                score3 = matrix[i - 1][j] + int(gap)
            else:
                score2 = matrix[i][j - 1]
                score3 = matrix[i-1][j]
            if (alignment == 'L'):
                matrix[i][j] = max(score1, score2, score3, 0) #Does not allow for negative numbers in local alignment
            else:
                matrix[i][j] = max(score1, score2, score3)

    return matrix
'''
buildDirectional- Creates the directional string for following the matrix to 
the best solution
'''
def buildDirectional(matrix, N, M, gapScore,alignment):
    directionalString = ''
    currentRow = N
    currentCol = M
    while(currentRow != 0 or currentCol != 0):

        if(alignment == 'L' and matrix[currentRow][currentCol] == 0): #returns dstring if local alignment is complete
            return directionalString
        if(currentRow == 0):
            directionalString = directionalString + ('H')
            currentCol = currentCol - 1
        elif(currentCol == 0):
            directionalString = directionalString + ('V')
            currentRow = currentRow - 1
        elif(matrix[currentRow][currentCol] == matrix[currentRow-1][currentCol] and alignment == 'S' and currentCol == M):
            directionalString = directionalString + ('V') #Moves semi local alignment away from terminal gaps
            currentRow = currentRow-1
        elif(matrix[currentRow][currentCol] == matrix[currentRow][currentCol-1] and alignment == 'S' and currentRow == N):
            directionalString = directionalString + ('H')
            currentCol = currentCol -1
        elif(matrix[currentRow][currentCol - 1] + int(gapScore) == matrix[currentRow][currentCol]):
            directionalString = directionalString + ('H')
            currentCol = currentCol - 1
        elif(matrix[currentRow - 1][currentCol] + int(gapScore) == matrix[currentRow][currentCol]):
            directionalString = directionalString + ('V')
            currentRow = currentRow - 1
        else:
            directionalString = directionalString + ('D')
            currentRow = currentRow - 1
            currentCol = currentCol - 1
    return directionalString
def buildAlignment(sequenceOne, sequenceTwo, directionalString):
    
    seq1Pos = len(sequenceOne)-1
    seq2Pos = len(sequenceTwo)-1
    dirPos = 0
    alignSeq1 = ''
    alignSeq2 = ''
    while(dirPos < len(directionalString)):
        
        if(directionalString[dirPos] == "D"): #Align sequence for match
            alignSeq1 = sequenceOne[seq1Pos] + alignSeq1
            alignSeq2 = sequenceTwo[seq2Pos] + alignSeq2
            seq1Pos = seq1Pos - 1
            seq2Pos = seq2Pos - 1
        elif(directionalString[dirPos] == "V"):#Place a space for a gap
            alignSeq1 = sequenceOne[seq1Pos] + alignSeq1
            alignSeq2 = '-' + alignSeq2
            seq1Pos = seq1Pos - 1
        else:
            alignSeq1 = '-' + alignSeq1 #Place a space for a gap
            alignSeq2 = sequenceTwo[seq2Pos] + alignSeq2
            seq2Pos = seq2Pos -1
        dirPos = dirPos + 1
    alSeq1List = list(alignSeq1)
    alSeq2List = list(alignSeq2)

    connectors, alignmentScore = createConnections(alSeq1List, alSeq2List)

    return alignSeq1, alignSeq2, connectors, alignmentScore
def createConnections(alSeq1List, alSeq2List):

    connectors = ""
    alignmentScore = 0
    for char1, char2 in zip(alSeq1List, alSeq2List):
        if char1 == char2:
            connectors = connectors + "|" #Create list for relationship between sequences
            alignmentScore = alignmentScore + 1

        else:
            if char1 == "-" or char2 == "-":
                connectors = connectors + " "
            else:
                connectors = connectors + "."

    return connectors, alignmentScore
def printAlignment(alignSeq1, alignSeq2, connectors, outfile):
    line = 0
    nucleotideCount = 1
    charPerLine = 50
    for line in range(line,len(alignSeq1),charPerLine):
        print(str(nucleotideCount) + " " + alignSeq1[line:line + charPerLine])
        outfile.write(str(nucleotideCount) + " " + alignSeq1[line:line + charPerLine]+ '\n')
        print(len(str(nucleotideCount)) * " " + " " + connectors[line:line + charPerLine])
        outfile.write(len(str(nucleotideCount)) * " " + " " + connectors[line:line + charPerLine]+ '\n')
        print(str(nucleotideCount) + " " + alignSeq2[line:line+charPerLine])
        outfile.write(str(nucleotideCount) + " " + alignSeq2[line:line+charPerLine])
        nucleotideCount = nucleotideCount + charPerLine
    return
def printStats(alignmentScore, SequenceOneLength, SequenceTwoLength, alignSeq1, alignSeq2, matrix, outfile):

    print("Match percentage is: " + str(alignmentScore/float(len(alignSeq1))*100))
    outfile.write("\nMatch percentage is: " + str(alignmentScore/float(len(alignSeq1))*100) + '\n')
    print("The alignment score is: " + str(matrix[SequenceOneLength][SequenceTwoLength]))
    outfile.write("The alignment score is: " + str(matrix[SequenceOneLength][SequenceTwoLength])+ '\n' )
    return
def main ():

    # Opens files
    infile1 = open('file1.txt', 'r')
    infile2 = open('file2.txt', 'r')
    outfile = open('ch3sk3out.txt', 'w')
    sequenceOne = readInFile(infile1)
    sequenceTwo = readInFile(infile2)

    sequenceOne,sequenceTwo = checkSequenceLength(sequenceOne,sequenceTwo)


    #Initialize matrix and prompt use for score preferences
    connectors = ''
    matrix = ''
    matrix = list(matrix)
    matrix = [[0 for col in range(len(sequenceTwo)+1)] for row in range(len(sequenceOne)+1)]
    
    print("This program aligns DNA sequences with global and semi global options.")
    gap = raw_input("Please enter gap score: ")
    misMatch = raw_input("Please enter mismatch score: ")
    match = raw_input("Please enter match score: ")
    alignment = raw_input("[G]lobal or [S]emi global alignment or [L]ocal: ")

    alignment = alignment.upper()

    matrix = buildMatrix(matrix,alignment,gap,misMatch,match,sequenceOne,sequenceTwo)
    
    #step 2 create directional string

    directionalString = ''
    maxValue = ''
    if(alignment == 'G' or alignment == 'S'):
        directionalString = buildDirectional(matrix, len(sequenceOne), len(sequenceTwo), gap,alignment)
    else:   
        maxValue = max([max(row) for row in matrix]) #Find the best alignment in matrix
        for rowIdx, rowValue in enumerate(matrix):
            if max(rowValue) == maxValue:
                maxRowIdx = rowIdx
                maxColIdx = rowValue.index(max(rowValue))
                directionalString = buildDirectional(matrix, maxRowIdx, maxColIdx, gap,alignment)
    #step 3 Build alignment using directional strings

    matrixOut = ''
    alignSeq1, alignSeq2, connectors, alignmentScore = buildAlignment(sequenceOne, sequenceTwo, directionalString)

    printAlignment(alignSeq1, alignSeq2, connectors, outfile)

    printStats(alignmentScore, len(sequenceOne), len(sequenceTwo), alignSeq1, alignSeq2, matrix, outfile)
    matrixOut = raw_input("Would you like to print the matrix? (1 yes, 2 no) ")
    if (matrixOut == "1"):
    	printMatrix(matrix,outfile)
    dStringOut = raw_input("Would you like to print the path string? (1 yes, 2 no) ")
    if (dStringOut == "1"):
        print (directionalString)
       	outfile.write('\n' + directionalString)
main()
