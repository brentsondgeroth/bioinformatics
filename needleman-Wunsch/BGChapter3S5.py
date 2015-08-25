'''
Name: Brent Gaither and Bill Moers
Description: This program uses the needleman wunch algorithm to align 
sequences globally or semi globally. 
bgwmCh3S5

Due Date: 1/22/15
'''

# Function that takes the lists that make up the matrix and format
# them to print in a table.
def printMatrix(matrix,outfile):
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print ('\n'.join(table))
    outfile.write('\n'.join(table))
# Function that builds a string that provides direction for tracing
# a path back through the matrix to build the best match.
def buildDirectional(matrix, N, M, gap):
    dString = ''
    currentRow = N
    currentCol = M
    gapScore = gap
    while(currentRow != 0 or currentCol != 0):
        if(currentRow == 0):
            dString = dString + ('H')
            currentCol = currentCol - 1
        elif(currentCol == 0):
            dString = dString + ('V')
            currentRow = currentRow - 1
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

    # Prompt user for input files
    '''
    data1 = raw_input("Enter name of wild-type file: ") 
    data2 = raw_input("Enter name of allele file: ")
    '''
    # Opens files
    infile1 = open('local1.txt', 'r')
    infile2 = open('local2.txt', 'r')
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

    #Intialize matrix and prompt use for score preferences
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
    alignment = input("[G]lobal or [S]emi global alignment: ")
    alignment = alignment.upper()
    
    #Step 1 build matrix
    matrix [0][0] = 0
    for i in range(1,N+1):
        if alignment == 'G':
            matrix [i][0] = matrix [i - 1][0] + int(gap)
        else:
            matrix [i][0] = 0
    for j in range(1,M +1 ):
        if alignment == 'G':
            matrix [0][j] = matrix[0][j - 1] + int(gap)
        else:
            matrix [0][j] = 0
    for i in range(1,N +1 ):
        for j in range(1,M+1):
            if (s1[i-1] == s2[j - 1]):
                score1 = matrix[i - 1][j - 1] + int(match)               
            else:
                score1 = matrix[i - 1][j - 1] + int(misMatch)
            
            score2 = matrix[i][j - 1] + int(gap)
            score3 = matrix[i - 1][j] + int(gap)
            
            matrix[i][j] = max(score1, score2, score3)

    #step 2 create directional string
    dString = ''
    alignment == 'G'
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
    lineCount = 1
    charPerLine = 50
    
    # Prints the alignment with a limit of 50 characters per line.
    for line in range(line,len(alSeq1),charPerLine):
        print(str(lineCount) + " " + alSeq1[line:line + charPerLine])
        outfile.write(str(lineCount) + " " + alSeq1[line:line + charPerLine]+ '\n')
        print(len(str(lineCount)) * " " + " " + connectors[line:line + charPerLine])
        outfile.write(len(str(lineCount)) * " " + " " + connectors[line:line + charPerLine]+ '\n')
        print(str(lineCount) + " " + alSeq2[line:line+charPerLine])
        outfile.write(str(lineCount) + " " + alSeq2[line:line+charPerLine])
        lineCount = lineCount + 1
    
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
