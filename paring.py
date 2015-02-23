
#how long the words will be, probably should remain the same as alignment
wordLength = 5

#Lower bound for inclusion in high scoring words, lower than alignment to allow for variance?
highScoreThreshold = 3

#Distance that two words must be under to group into a segment pair
segPairWordMaxDist = 12

#How far around the high score segment should be included?
paringBufferLength = 2

#Takes a sequence (or mean) and returns the list of words of length wordLength
#Used lots, copied here for access to wordLength
def wordify (seq):
	seqWords = []
	for i in range((len(seq)-wordLength)):
		seqWords += [seq[i:(i+wordLength)]]
	return seqWords

def scoreMeanWord(meanWord):
	score = 0
	for probArray in meanWord:
		score += max(probArray)
	return score

#Perhaps account for small gaps in a single motif here?
#Returns means from the paring data
def processSegments(segments, mean):
	print segments
	newMeans = []
	for (score,i,j) in segments:
		i1 = max(0, i - paringBufferLength)
		i2 = min(len(mean),j+wordLength+paringBufferLength)
		newMeans += [mean[i1:i2]]
	# print str(newMeans) + '\n'
	return newMeans

#Somewhat similar to alignment, generates high scoring fragments of centroids and branches them off as means
def pare(mean):
	meanWords = wordify(mean)
	numMeanWords = len(meanWords)
	scores = []
	for m in range(numMeanWords):
		score = scoreMeanWord(meanWords[m])
		scores += [score]
	segments = []
	i = 0
	while i < numMeanWords:
		if scores[i] > highScoreThreshold:
			currSegScore = 0
			j = 0
			#FIX THIS - what if there is a non-high score word in the middle of a motif?
			while j < segPairWordMaxDist and i+j < numMeanWords and scores[i+j] > highScoreThreshold:
				currSegScore += scores[i+j]
				j += 1
			#So stores it as score, beginning word, ending word
			#Score currently isn't used, but may be
			segments += [(currSegScore,i,i+j)]
			i += j
		i += 1
	#Segments need to turn into means that account for buffer distances and the length of meanWords
	return processSegments(segments, mean)

#Try paring without meanWords?
def pare2(mean):
	return


#Currently the algorithm pares every mean every iteration, maybe too much?
def paredMeans(means):
	newMeans = []
	for mean in means:
		newMeans += pare(mean)
	return newMeans