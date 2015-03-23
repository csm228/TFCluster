
import math

#how long the words will be, probably should remain the same as alignment
wordLength = 5

#Lower bound for inclusion in high scoring words, lower than alignment to allow for variance?
highScoreThreshold = 3

#Distance that two words must be under to group into a segment pair
segPairWordMaxDist = 12

#How far around the high score segment should be included?
paringBufferLength = 2

#How large of a standard deviation in mean alignment length should cause paring?
paringDevAllowance = 1.5

#How much difference from the mode length should be removed?
paringLengthSimilarity = 1

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
def originalPare(mean):
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



#Returns a mean from a selected alignment datum
def processAlignment((i,m,score,length), mean):
	print (i,m,score,length)
	i1 = max(0, m - paringBufferLength)
	i2 = min(len(mean),m+length+paringBufferLength)
	return mean[i1:i2]

#Used iteratively in the paring function to guess new means,
#Essentially just a variance function
def lengthVar(lengths):
	numLengths = float(len(lengths))
	sumLengths = 0
	for length in lengths:
		sumLengths += length
		lengths += [length]
	avgLength = sumLengths/float(numPeaks)
	sumVar = 0
	for (i,m,score,length) in lengthAlignments:
		sumVar += (avgLength + float(length))**2
	lengthVariance = sumVar / float(numPeaks)
	return lengthVariance

#Somewhat similar to alignment, generates high scoring fragments of centroids and branches them off as means
#Now uses mean alignment length to guess length and position of new means
#For a combined length and conserved alignment, 
def pare((mean, targetLength), lengthAlignments):
	#Odd artifact here, the targetLength is from alignments on the previous clustering,
	# so they don't match the new score. Take out the extra alignmentMatrix generation?
	#FIXED with summing new memoized alignments, but may want to check this
	numPeaks = len(cluster)-1
	if numPeaks <= 0:
		return []
	lengths = []
	for (i,m,score,length) in lengthAlignments:
		lengths += [length]
	stdDev = math.sqrt(lengthVar(lengths))
	newMeans = []
	#Currently this DOESN'T EVEN CHECK TARGETLENGTH
	while stdDev >= paringDevAllowance:
		#finds the most common length, the mode
		modeLength = max(set(lengths), key=lengths.count)
		selectedAlignment = (0,0,0,0)
		while k < len(lengthAlignments):
			#length = lengths[i]
			(i,m,score,length) = lengthAlignments[i]
			if length == modeLength:
				selectedAlignment = lengthAlignments.pop(i)
				lengths.pop(i)
			else:
				k += 1
		newMeans += [processAlignment(selectedAlignment,mean)]
		if len(lengthAlignments) <= 0:
			return newMeans
		stdDev = math.sqrt(lengthVar(lengths))
	return newMeans

#Currently the algorithm pares every mean every iteration, maybe too much?
def paredMeans(means,varAlignmentMatrix):
	newMeans = []
	for j in range(len(means)):
		newMeans += pare(means[j],alignmentMatrix[j])
	return newMeans

# def paredMeans(means,alignmentMatrix):
# 	newMeans = []
# 	for mean in means:
# 		newMeans += originalPare(mean)
# 	return newMeans