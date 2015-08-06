
import bioToolsP


pickAgain = True

print("This program has 5 bioinformatics tools to use. The needleman wunch tool aligns " +
	"DNA sequences globally or semi globally or locally. The gene Probablity tool gives "
	+ "the probability that a eukaryotic or prokaryotic sequence is a gene. The substitution matrix "+
	"tool creates a matrix that gives the probablity of nucleotides aligning with each other in a substitution." + 
	"The hidden markov model predicts the probability of a eukaryotic sequence. The phylogenetic tree "+
	"creates the most likely tree and outputs the tree in Newick format.\n For all files include suffix when "+
	"inputting file names.") 

while (pickAgain):

	userChoice = raw_input("Select a tool to use \n 1. Needleman Wunch Algorithm \n 2. Gene Probablity (based on user " +
		"input for promoter sequences\n 3. Substitution Matrix \n 4. Hidden Markov Model(based on a predefined model)\n "+
		"5. Phylogenetic Tree\n 6. Exit\n")

	if userChoice == '1':
		bioToolsP.needleman()
	elif userChoice == '2':
		bioToolsP.geneProbablity()
	elif userChoice == '3':
		bioToolsP.substitutionMatrix()
	elif userChoice == '4':
		bioToolsP.markovModel()
	elif userChoice == '5':
		bioToolsP.phylogeneticTree()
	else:
		pickAgain = False

