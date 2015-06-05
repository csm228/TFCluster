# 
# USE: python motifMaker.py <number of seeds> <number of footprints> <min footprint length> <max footprint length> <outfile name>
# OUTPUT: Footprint sequence, sequence position(number), motif start position, Seed info (list)
import random
import sys


#how long should the generated motifs be? (Change to an input for a command/bash?)
maxMotifLength = 7
# motifLength = int(sys.argv[1])

#alternative mass generation
minMotifLength = 5
maxMotifLength = 9

bases = ['A','T','G','C']

def recMotifDictBuilder(listMotifs,currMotif,lengthCount):
	if lengthCount <= 0:
		listMotifs = listMotifs + currMotif + "\": 0, \""
		return listMotifs
	else:
		for base in bases:
			listMotifs = recMotifDictBuilder(listMotifs, currMotif + base, lengthCount - 1)
		return listMotifs


# main. Get info, generate motif tuple, write footprints to output.
# motifLength = int(sys.argv[1])
outfile = open('variableSet.py', 'w')
for i in range(minMotifLength,maxMotifLength):
	outfile.write("\n\nmotifDictL" + str(i) + " = ")
	motifDict = recMotifDictBuilder("","",i)
	motifDict = "{\"" + motifDict[:-3] + "}"
	outfile.write(motifDict)
outfile.close()