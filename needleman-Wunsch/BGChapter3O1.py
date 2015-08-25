'''
Name: Brent Gaither 
Description: This program uses the needleman wunch algorithm to align 
sequences globally or semi globally or locally 
bgwmCh3O1
Due Date: 1/27/15
'''

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
def buildMatrix(matrix,alignment,gap,misMatch,match,N,M,s1,s2):

    #Step 1 build matrix
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
            if (s1[i-1] == s2[j - 1]):
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
    dString = ''
    currentRow = N
    currentCol = M
    while(currentRow != 0 or currentCol != 0):

        if(alignment == 'L' and matrix[currentRow][currentCol] == 0): #returns dstring if local alignment is complete
            return dString
        if(currentRow == 0):
            dString = dString + ('H')
            currentCol = currentCol - 1
        elif(currentCol == 0):
            dString = dString + ('V')
            currentRow = currentRow - 1
        elif(matrix[currentRow][currentCol] == matrix[currentRow-1][currentCol] and alignment == 'S' and currentCol == M):
            dString = dString + ('V') #Moves semi local alignment away from terminal gaps
            currentRow = currentRow-1
        elif(matrix[currentRow][currentCol] == matrix[currentRow][currentCol-1] and alignment == 'S' and currentRow == N):
            dString = dString + ('H')
            currentCol = currentCol -1
        elif(matrix[currentRow][currentCol - 1] + int(gapScore) == matrix[currentRow][currentCol]):
            dString = dString + ('H')
            currentCol = currentCol - 1
        elif(matrix[currentRow - 1][currentCol] + int(gapScore) == matrix[currentRow][currentCol]):
            dString = dString + ('V')
            currentRow = currentRow - 1
        else:
            dString = dString + ('D')
            currentRow = currentRow - 1
            currentCol = currentCol - 1
    return dString

def main ():

    # Opens files
    infile1 = open('ch5sequence1.txt', 'r')
    infile2 = open('ch5sequence2.txt', 'r')
    outfile = open('ch3sk3out.txt', 'w')
    
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
    temp = ""
    if (len(s1) > len(s2)): #ensure long sequence is s2
        temp = s1
        s1 = s2
        s2 = temp

    #Initialize matrix and prompt use for score preferences
    N = len(s1)
    M = len(s2)
    connectors = ''
    matrix = ''
    matrix = list(matrix)
    matrix = [[0 for col in range(M+1)] for row in range(N+1)]
    print("This program aligns DNA sequences with global and semi global options.")
    gap = input("Please enter gap score: ")
    misMatch = input("Please enter mismatch score: ")
    match = input("Please enter match score: ")
    alignment = input("[G]lobal or [S]emi global alignment or [L]ocal: ")
    alignment = alignment.upper()

    buildMatrix(matrix,alignment,gap,misMatch,match,N,M,s1,s2)
    #step 2 create directional string
    dString = ''
    maxValue = ''
    if(alignment == 'G' or alignment == 'S'):
        dString = buildDirectional(matrix, N, M, gap,alignment)
    else:   
        maxValue = max([max(row) for row in matrix]) #Find the best alignment in matrix
        for rowIdx, rowValue in enumerate(matrix):
            if max(rowValue) == maxValue:
                maxRowIdx = rowIdx
                maxColIdx = rowValue.index(max(rowValue))
                N = maxRowIdx
                M = maxColIdx
                dString = buildDirectional(matrix, N, M, gap,alignment)
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
    nucleotideCount = 1
    charPerLine = 50
    
    # Prints the alignment with a limit of 50 characters per line.
    for line in range(line,len(alSeq1),charPerLine):
        print(str(nucleotideCount) + " " + alSeq1[line:line + charPerLine])
        outfile.write(str(nucleotideCount) + " " + alSeq1[line:line + charPerLine]+ '\n')
        print(len(str(nucleotideCount)) * " " + " " + connectors[line:line + charPerLine])
        outfile.write(len(str(nucleotideCount)) * " " + " " + connectors[line:line + charPerLine]+ '\n')
        print(str(nucleotideCount) + " " + alSeq2[line:line+charPerLine])
        outfile.write(str(nucleotideCount) + " " + alSeq2[line:line+charPerLine])
        nucleotideCount = nucleotideCount + charPerLine
    
    print("Match percentage is: " + str(alignmentScore/float(len(alSeq1))*100))
    outfile.write("\nMatch percentage is: " + str(alignmentScore/float(len(alSeq1))*100) + '\n')
    print("The alignment score is: " + str(matrix[N][M]))
    outfile.write("The alignment score is: " + str(matrix[N][M])+ '\n' )
    matrixOut = input("Would you like to print the matrix? (1 yes, 2 no) ")
    if (matrixOut == "1"):
    	printMatrix(matrix,outfile)
    dStringOut = input("Would you like to print the path string? (1 yes, 2 no) ")
    if (dStringOut == "1"):
        print (dString)
       	outfile.write('\n' + dString)
main()
