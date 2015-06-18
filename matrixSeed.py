# import align

import variableSet
# import random

#how long the sequences in the matrix will be
#May want to set an input function to recalculate parameters in the event that the length changes?
wordLength = 6
motifDict = variableSet.motifDictL6

#The additional blank (0.25 each) character positions added to either side of new seeds
bufferLength = 1
blankProb = [0.25,0.25,0.25,0.25]

#Lower bound for inclusion in high scoring words
highScoreThreshold = 3.5




# The matrix array of characters is a list of [probA;probT;probG;probC] probLists
def initProb (character):
	if character == 'A':
		return [1.0,0.0,0.0,0.0]
	elif character == 'T':
		return [0.0,1.0,0.0,0.0]
	elif character == 'G':
		return [0.0,0.0,1.0,0.0]
	elif character == 'C':
		return [0.0,0.0,0.0,1.0]
	else:
		print "Incorrect string in peaks, implement error handling"

#abstracts a sequence seed into a mean
#CAREFUL, currently the only thing with action adding mean length :P
def abstractSeed (seedSeq):
	matrix = []
	matrix += [blankProb] * bufferLength
	for character in seedSeq:
		matrix += [initProb(character)]
	matrix += [blankProb] * bufferLength
	#if the words are high frequency, we can already assume they are means
	return (matrix, wordLength, 0)

#Takes a sequence (or mean) and returns the list of words of length wordLength
#Used in the alignment function
def wordify (seq):
	seqWords = []
	for i in range((len(seq)-wordLength+1)):
		seqWords += [seq[i:(i+wordLength)]]
	return seqWords

	
#IDEA: A score-based (not purely frequency based) version of matrix membership? (follows)

#The matrix array of characters is a list of [probA;probT;probG;probC] lists
def compare (character, probArray):
	if character == 'A':
		return probArray[0]
	elif character == 'T':
		return probArray[1]
	elif character == 'G':
		return probArray[2]
	elif character == 'C':
		return probArray[3]
	else:
		print "Incorrect string in peaks, implement error handling"

#assume both have the same length, generally = wordLength
def scorePair (seqWord,meanWord):
	#if misalignments become an issue, use weighted mean?
	score = 0
	for i in range(len(seqWord)):
		if seqWord[i] == meanWord[i]:
			score += 1
	return score


def scoreSeedMatrix(peaks):
	for peak in peaks:
		seq = peak[0]
		seqWords = wordify(seq)
		for word in seqWords:
			for key in motifDict:
				score = scorePair(key,word)
				if score > highScoreThreshold:
					motifDict[key] = motifDict[key] + score
	return

#Now for purely seed and calculation methods: perhaps base the number of seeds on the matrix results?
#Other ideas for seed number (k) changing?

# def tally(keySeq):
# 	variableSet.motifDict[keySeq] = variableSet.motifDict[keySeq] + 1

def tallySeedMatrix(peaks):
	for peak in peaks:
		seq = peak[0]
		seqWords = wordify(seq)
		for keySeq in seqWords:
			motifDict[keySeq] = motifDict[keySeq] + 1
	return

#Assumes numSeeds <= length of the dictionary. A perfectly reasonable assumption.
def firstNSeeds(numSeeds):
	seeds = []
	seqFreq = motifDict.items()
	print seqFreq[0]
	seqFreq.sort(key = lambda seg: seg[1], reverse=True)
	print seqFreq[0]
	for i in range(numSeeds):
		(key,value) = seqFreq[i]
		seeds += [abstractSeed(key)]
	return seeds

#Assumes that the initial mean guess won't be larger that the number of peaks
#May implement so that it is the only function to use subsampling,
# assumes the initial k guess will reduce the number of sequences in any cluster significantly

def pickInitSeeds (peaks, numMeans):
	scoreSeedMatrix(peaks)
	return firstNSeeds(numMeans)

# def pickInitSeeds (peaks, numMeans):
# 	tallySeedMatrix(peaks)
# 	return firstNSeeds(numMeans)