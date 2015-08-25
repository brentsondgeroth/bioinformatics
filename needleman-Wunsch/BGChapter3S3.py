'''
Name: Brent Gaither and Bill Moers
Description: This program uses the needleman wunch algorithm 
Due Date: 1/20/15
'''

import textwrap

# Takes in a sequence and makes the text wrap to a limit
# of 20 characters.
def wrapText(seq):
    width = 20
    list1 = textwrap.wrap(seq, width)
    for element in list1:
        print (element)
    return

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

    # Prompt user for input
    '''
    data1 = raw_input("Enter name of wild-type file: ") 
    data2 = raw_input("Enter name of allele file: ")
    '''
    # Opens files
    infile1 = open('file3.txt', 'r')
    infile2 = open('file4.txt', 'r')
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

    N = len(s1)
    M = len(s2)
    connectors = ''
    matrix = ''
    matrix = list(matrix)
    matrix = [[0 for col in range(M+1)] for row in range(N+1)]
    gap = input("Please enter gap score: ")
    misMatch = input("Please enter mismatch score: ")
    match = input("Please enter match score: ")

    #Step 1 build matrix
    matrix [0][0] = 0
    for i in range(1,N - 1):
        matrix [i][0] = matrix [i - 1][0] + int(gap)
    for j in range(1,M - 1):
        matrix [0][j] = matrix[0][j - 1] + int(gap)
    for i in range(1,N - 1):
        for j in range(1,M - 1):
            if (s1[i-1] == s2[j - 1]):
                score1 = matrix[i - 1][j - 1] + int(match)               
            else:
                score1 = matrix[i - 1][j - 1] + int(misMatch)
            
            score2 = matrix[i][j - 1] + int(gap)
            score3 = matrix[i - 1][j] + int(gap)
            
            matrix[i][j] = max(score1, score2, score3)
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
        if(dString[dirPos] == "D"):
            alSeq1 = s1[seq1Pos] + alSeq1
            alSeq2 = s2[seq2Pos] + alSeq2
            seq1Pos = seq1Pos - 1
            seq2Pos = seq2Pos - 1
        elif(dString[dirPos] == "V"):
            alSeq1 = s1[seq1Pos] + alSeq1
            alSeq2 = '-' + alSeq2
            seq1Pos = seq1Pos - 1
        else:
            alSeq1 = '-' + alSeq1
            alSeq2 = s2[seq2Pos] + alSeq2
            seq2Pos = seq2Pos -1

        dirPos = dirPos + 1

    alSeq1List = list(alSeq1)
    alSeq2List = list(alSeq2)

    for char1, char2 in zip(alSeq1List, alSeq2List):
        if char1 == char2:
            connectors = connectors + "|"
            alignmentScore = alignmentScore + 1

        else:
            if char1 == "-" or char2 == "-":
                connectors = connectors + " "
            else:
                connectors = connectors + "."
    
    line = 0
    lineCount = 1
    charPerLine = 20
    
    for line in range(line,len(alSeq1),charPerLine):
        print (str(lineCount) + " " + alSeq1[line:line + charPerLine])
        print(len(str(lineCount)) * " " + " " + connectors[line:line + charPerLine])
        print(str(lineCount) + " " + alSeq2[line:line+charPerLine])
        lineCount = lineCount + 1

    #print (alSeq1)
    #print (alSeq2)
main()


        while(True):
            if(dString[dirPos] == "D"):
                break
            if(dString[dirPos] == "V"):
                alSeq1 = s1[seq1Pos] + alSeq1
                alSeq2 = '-' + alSeq2
                seq1Pos = seq1Pos - 1
            else:
                alSeq1 = '-' + alSeq1
                alSeq2 = s2[seq2Pos] + alSeq2
                seq2Pos = seq2Pos -1
            dirPos = dirPos + 1
