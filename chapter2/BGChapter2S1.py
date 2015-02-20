'''
Names:Matthew Kachlik and Brent Gaither
Description: This program takes information from the user to proecess the 
two sequences of DNA then comapres the amino acid sequence looking for 
differences then writes to an outputfile
Due Date: 1/15/15
'''
import sys


'''
    transcription
    takes in the a sequence then replaces all T's with U's 
    then returns the new string 
'''
def transcription(seq):
    seq = seq.replace('T', 'U')
    return seq
'''
    comparision
    takes in two sequenes then compares them looking for differences
    if there are errors the position and error is written out to a file
'''
def comparision(seq1,seq2,typeOfComparison):
    outfile = open('ch2sk1out.txt', 'w')
    
    seq1 = seq1
    posErrors = ""
    #check for mutations 
    print (typeOfComparison)
    mctr = 0

    posErrors = "Unmatched " + typeOfComparison + " at "
    for iterator in range(0,len(seq1)):
        #when you find a mutation write to the file where the error is
        if(seq1[iterator] != seq2[iterator]):
             posErrors = posErrors + " " + str(iterator+1)
             mctr += 1
    if(mctr == 0):
        print ("no mismatches found - strings are identical")
    else:
        print(posErrors)
        print(seq1)
        print(seq2)
        outfile.write(seq1 + " " + str(iterator+1) + " " + seq2)
    outfile.close()
'''
    translation
    takes in two seqs and a hash table of codes then converts 
    the two sequences of genetic code into list then uses the hash
    table to make an amino acid sequence then calls comparision 
    with both amino acid sequences
'''
def translation(seq1,seq2,codes):
    seq1List = list(seq1)
    seq2List = list(seq2)
    aminoSeq1 = ""
    aminoSeq2 = ""
    # Step 3: Translation do for both
    for iterator in range(0,len(seq1List)-1,3):
        key1 = ''.join(seq1List[iterator:iterator+3])
        aminoSeq1 = aminoSeq1 +codes[key1]
        key2 = ''.join(seq2List[iterator:iterator+3])
        aminoSeq2 = aminoSeq2 +codes[key2]

    #once translated compare them 
    comparision(aminoSeq1,aminoSeq2,"*Amino Acid sequence*")
'''
    complimentSequence
    takes in a seq then finds the compliment by 
    replacing all the A's with T's G's with C's T's 
    with A's with C's with G's then reutrns the sequence
'''
def complimentSequence(seq):
    #using temp variables to change T's for A's and G's for C's
    seq = seq.replace("A","%temp%").replace("T","A").replace("%temp%","T")
    seq = seq.replace("C","%temp%").replace("G","C").replace("%temp%","G")
    return seq

def main():

    # Step 1: Initialization and Read in Sequences
    codes = dict(GCA='A', GCC='A', GCG='A',GCU='A',
             AGA='R', AGG='R', CGA='R', CGC='R', CGG='R', CGU='R',
             AAU='N', AAC='N',
             GAC='D', GAU='D',
             UGC='C', UGU='C',
             GAA='E', GAG='E',
             CAA='Q', CAG='Q',
             GGA='G', GGC='G', GGG='G', GGU='G',
             CAC='H', CAU='H',
             AUA='I', AUC='I', AUU='I',
             CUA='L', CUC='L', CUG='L', CUU= 'L', UUA='L', UUG='L',
             AAA='K',AAG='K', AUG='M',
             UUC='F', UUU='F',
             CCA='P', CCC='P', CCG='P', CCU='P', 
             AGC='S', AGU='S', UCA='S', UCC='S', UCG= 'S', UCU='S', 
             ACA='T', ACC='T', ACG='T', ACU='T',
             UGG='M',
             UAC='Y', UAU='Y', GUA='V', GUC='V', GUG='V', GUU='V',
             UAA='*', UAG='*', UGA='*')

    infile1 = open('wildtype.txt', 'r')
    infile2 = open('mutantA.txt', 'r')
    
    #step 2 read fasta formated sequences
    seq1 = ""
    infile1.readline()  # bypass > header line
    for line in infile1:
        line = line.replace('\n', '')
        line = line.upper()
        seq1 = seq1 + line

    seq2 = ""
    infile2.readline()  # bypass > header line
    for line in infile2:
        line = line.replace('\n', '')
        line = line.upper()
       
        seq2 = seq2 + line
    #quit if they arent the same length
    if len(seq1) != len(seq2):
        print ("Sequences invalid - unequal length")
        sys.exit()

    # step 3 get format of the sequences and process sequencese to
    # correct format
    print ("This program compares DNA squences to find mutations")
    print ("Enter what type of template the sequences are")
    isTemplate = input("1) For non template sequences\n2) For template sequences\n")
    if isTemplate == "2":
        seq1 = complimentSequence(seq1)
        seq2 = complimentSequence(seq2)
    print ("Enter the direction of the strand")
    direction = input("1) 5' to 3' \n2) 3' to 5'\n")
    if direction == "2":
        seq1 = seq1[::-1]
        seq2 = seq2[::-1]
    print ("Enter The type of comparison you would like to make")

    user_input = input("1) For Direct DNA comparison\n2) For Translation then a comparision\n3) For Both comparisions\n")
    # step 4 do the comparision of the sequences
    '''
        This if statement call methods for more preprocessing
        of the sequences if neccasary 
    '''
    if user_input == "1":
        comparision(seq1,seq2,"*DNA Sequence*")
    elif user_input == "2":
        seq1 = transcription(seq1)
        seq2 = transcription(seq2)
        translation(seq1,seq2,codes)
    elif user_input == "3":
        comparision(seq1,seq2,"*DNA Squence*")
        seq1 = transcription(seq1)
        seq2 = transcription(seq2)
        translation(seq1,seq2,codes)
    #step 5 close up the files before exiting
    infile1.close()
    infile2.close()
main()
