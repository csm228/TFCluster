
import math

#how long the words will be, probably should remain the same as alignment
wordLength = 5

#Lower bound for inclusion in high scoring words, lower than alignment to allow for variance?
highScoreThreshold = 3

#Distance that two words must be under to group into a segment pair
segPairWordMaxDist = 13

#How far around the high score segment should be included?
paringBufferLength = 1
blankProb = [0.25,0.25,0.25,0.25]

#How large of a standard deviation in mean alignment length should cause paring?
paringDevAllowance = 1.0

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
#Returns means from the paring data, modified for targetLength addition
def processSegments(segments, mean):
	print segments
	newMeans = []
	for (score,i,j) in segments:
		overflowL = paringBufferLength - i
		overflowR = j + paringBufferLength - len(mean)
		newMean = []
		if overflowL > 0:
			if overflowR > 0:
				newMean = ([blankProb] * overflowL) + mean + ([blankProb] * overflowR)
			else:
				newMean = ([blankProb] * overflowL) + mean[:j]
		elif overflowR > 0:
			newMean = mean[i:] + ([blankProb] * overflowR)
		else:
			newMean = mean[i:j]
		newMeans += [(newMean,-1)]
	# print str(newMeans) + '\n'
	return newMeans

#Somewhat similar to alignment, generates high scoring fragments of centroids and branches them off as means
#Modified for use with (sequence, score) means - unboxing and return type
def originalPare((mean,targetLength)):
	meanWords = wordify(mean)
	numMeanWords = len(meanWords)
	scores = []
	for m in range(numMeanWords):
		score = scoreMeanWord(meanWords[m])
		scores += [score]
	segments = []
	i = 0
	while i < numMeanWords:
		currSegScore = scores[i]
		if currSegScore > highScoreThreshold:
			j = 1  # want to start at the next word
			#FIX THIS - what if there is a non-high score word in the middle of a motif? (actually, that's pretty unlikely)
			while j < segPairWordMaxDist and i+j < numMeanWords and scores[i+j] > highScoreThreshold:
				currSegScore += max(mean[i+j])
				j += 1
			#So stores it as score, beginning word, ending word
			#Score currently isn't used, but may be
			segments += [(currSegScore,i,i+j+wordLength)]
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
	avgLength = sumLengths/float(numLengths)
	sumVar = 0
	for length in lengths:
		sumVar += (avgLength - float(length))**2
	lengthVariance = sumVar / float(numLengths)
	return lengthVariance


#BELOW: SERIOUS BUGCHECK MAY BE NECESSARY, as with originaPare()

#Somewhat similar to alignment, generates high scoring fragments of centroids and branches them off as means
#Now uses mean alignment length to guess length and position of new means
#For a combined length and conserved alignment, 
def pare((mean, targetLength), lengthAlignments):
	#Odd artifact here, the targetLength is from alignments on the previous clustering,
	# so they don't match the new score. Take out the extra alignmentMatrix generation?
	#FIXED with summing new memoized alignments, but may want to check this
	numPeaks = len(lengthAlignments)
	if numPeaks <= 0:
		return []
	lengths = []
	for (i,m,score,length) in lengthAlignments:
		lengths += [length]
	stdDev = math.sqrt(lengthVar(lengths))
	newMeans = []
	#Currently this DOESN'T check or utilize TARGETLENGTH
	if stdDev <= paringLengthSimilarity:
		modeLength = max(set(lengths), key=lengths.count)
		l = 0
		while 0 <= l < numPeaks:
			(i,m,score,length) = lengthAlignments[l]
			if length == modeLength:
				newMeans += [(processAlignment(lengthAlignments[l],mean),-1)]
				l = -1
			else:
				l += 1
		return newMeans
	while stdDev >= paringDevAllowance:
		#finds the most common length, the mode
		modeLength = max(set(lengths), key=lengths.count)
		#Temp value, selectedAlignment should always be changed, as something has the mode length
		selectedAlignment = (0,0,0,0)
		k = 0
		while k < len(lengthAlignments):
			#length = lengths[i]
			(i,m,score,length) = lengthAlignments[k]
			if length == modeLength:
				#removes the values for a second iteration
				selectedAlignment = lengthAlignments.pop(k)
				lengths.pop(k)
			elif length - paringLengthSimilarity <= length <= length + paringLengthSimilarity:
				#removes close lengths too
				lengthAlignments.pop(k)
				lengths.pop(k)
			else:
				k += 1
		newMeans += [(processAlignment(selectedAlignment,mean),-1)]
		if len(lengthAlignments) <= 0:
			return newMeans
		stdDev = math.sqrt(lengthVar(lengths))
	return newMeans

# #Currently the algorithm pares every mean every iteration, maybe too much?
# def paredMeans(means,varAlignmentMatrix):
# 	newMeans = []
# 	for j in range(len(means)):
# 		#make sure the pare function doesn't change the alignment matrix
# 		newMeans += pare(means[j],list(varAlignmentMatrix[j]))
# 	return newMeans

#Uses originalPare, but unboxes 
def paredMeans(means,alignmentMatrix):
	newMeans = []
	for mean in means:
		newMeans += originalPare(mean)
	return newMeans