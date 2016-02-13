'''
Name: Brent Gaither
Description: This program takes information from the user to process the 
two sequences of DNA then compares the sequences to find the best fit. 
This program also allows for sequences to be unequal length.
Due Date: 1/20/15
'''
import sys

def readInFile(infile):
    sequence = ""
    infile.readline()  # bypass > header line
    for line in infile:
        line = line.replace('\n', '')
        sequence = sequence + line
        sequence = sequence.upper()
    infile.close()
    return sequence


'''
    comparison
    takes in two sequences then compares them looking for differences
    if there are errors the position and error is written out to a file
'''
def comparison(seq1,seq2):
        
    counterComp = 0 #

    for iterator in range(0,len(seq1)):
        #find matching nucleotides and add total 
        if(seq1[iterator] == seq2[iterator]):
             counterComp += 1
    if(counterComp == 0):
        print ("no mismatches found - sequences are identical")
    else:        
        print(seq1)
        print(seq2)
        print("Sequences are identical length with " + str(counterComp) + " matches out of " 
            + str(len(seq1)) + " nucleotides")
    '''
    comparisonUnequal
    takes in two sequences then compares them looking for differences and 
    best fit for deletions. 
    '''
def comparisonUnequal(seq1,seq2):
    difference = 0
    deletion =''

    bestFit = list('')
    counter = 0
    if(len(seq1) < len(seq2)):
        seq1, seq2 = switchSequences(seq1,seq2)

    difference = len(seq1) - len(seq2)
    shortSeq = len(seq2)

    #create spaces for deletion
    for i in range(0,difference):
        deletion = deletion + '-'   

    #iterate through sequence moving deletion down each nucleotide
    for iterator in range(0,shortSeq+1):
        counter = 0 #resets counter to 0
        newSeq2 = seq2[0:iterator] + deletion + seq2[iterator:shortSeq] #move in deletion into iterator position
        #count number of matching nucleotides put in list
        for i in range(0,shortSeq+1):
            if seq1[i] == newSeq2[i]:
                counter = int(counter)+ 1
        bestFit.append(int(counter)) #add number of matching nucleotides to list
    print("The best matching sequence is...")
    print(seq1)
    print(seq2[0:bestFit.index(max(bestFit))] + deletion + seq2[bestFit.index(max(bestFit)):shortSeq])
    print("There are " + str(max(bestFit)) + " matching nucleotides")
    print("A deletion of " + str(difference) + " nucleotide(s) occurred at nucleotide(s) " 
        + str(bestFit.index(max(bestFit))+1) + "-"+ str(bestFit.index(max(bestFit))+difference))

def switchSequences(seq1,seq2):
    temp = seq1
    seq1 = seq2
    seq2 = temp
    return seq1, seq2

'''
    complimentSequence
    takes in a seq then finds the compliment by 
    replacing all the A's with T's G's with C's T's 
    with A's with C's with G's then returns the sequence
'''

def complimentSequence(seq):
    #using temp variables to change T's for A's and G's for C's
    seq = seq.replace("A","%temp%").replace("T","A").replace("%temp%","T")
    seq = seq.replace("C","%temp%").replace("G","C").replace("%temp%","G")
    return seq

def main():

    # Step 1: Initialization and Read in Sequences

    infile1 = open('mutantB copy.txt', 'r')
    infile2 = open('mutantB.txt', 'r')
    
    seq1 = readInFile(infile1)
    seq2 = readInFile(infile2)

    # step 2 get format of the sequences and process sequences to
    # correct format
    print ("This program compares DNA sequences to find mutations the first strand will " +
        "be considered the reference strand")
    print ("Enter what type of template the sequences are")
    isTemplate = raw_input("1) For non template sequences\n2) For template sequences\n")

    print ("Enter the direction of the strand")
    direction = raw_input("1) 5' to 3' \n2) 3' to 5'\n")

    if isTemplate == "2":
        seq1 = complimentSequence(seq1)
        seq2 = complimentSequence(seq2)

    if direction == "2":
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]
    
    # step 3 do the comparison of the sequences
    if(len(seq1) == len(seq2)):
        comparison(seq1,seq2)        
    else:
        comparisonUnequal(seq1,seq2)

    #step 4 close up the files before exiting
    infile1.close()
    infile2.close()
main()