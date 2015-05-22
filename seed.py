import align

import random

def calcSubSize (numMeans, numPeaks):
	#deal with the case where the subsample or the number of peaks left
	#is less than the number of new guessed clusters [+5 below]
	subSize = numPeaks//10 + 5
	return min(subSize,numPeaks) #This feels horrible - deal with a lack of outliers ELSEWHERE

#Takes a random peak and removes it from the list of peaks,
#returns (peak, remaining peaks) as a tuple
def sample (peaks):
	index = random.randrange(len(peaks))
	return (peaks.pop(index),peaks)

def minScore (peak, meanWordsList):
	#align gives a tuple of index, score, so take the score
	minVal = align.score_of_align(peak, meanWordsList[0])
	for meanWords in meanWordsList[1:]:
		minVal = min(minVal, align.score_of_align(peak,meanWords))
	return minVal

#Finds the distance to the closest mean for each peak,
#by comparing the distance to the newest mean with the value
#for closest mean distance memoized from the previous run
def newPeakDistances (mean,peaks,peakDistances):
	meanWords = align.wordify(mean)
	for i in range(len(peakDistances)):
		score = align.score_of_align(peaks[i],meanWords)
		peakDistances[i] = max(peakDistances[i],score)
	return peakDistances

#As an alternative, generate all possible sequences of a given length (7?), run k++ on the most common, then select representative sequences


def kPlusPlus (means, peaks, peakDistances):
	print "kPlusPlus"
	newestMean = means[len(means)-1]
	peakDistances = newPeakDistances(newestMean,peaks,peakDistances)
	seedIndex = 0
	minVal = peakDistances[0]
	for i in range(1,len(peaks)):
		#Bias? should this be randomized instead?
		tempScore = peakDistances[i]
		if tempScore < minVal:
			minVal = tempScore
			seedIndex = i
	seed = peaks.pop(seedIndex)
	# print seed
	peakDistances.pop(seedIndex)
	return (seed,peaks,peakDistances)

# def kPlusPlus (means, peaks):
# 	print "kPlusPlus"
# 	meanWordsList = []
# 	for mean in means:
# 		meanWordsList += [align.wordify(mean)]
# 	seedIndex = 0
# 	minVal = minScore(peaks[0], meanWordsList)
# 	for i in range(1,len(peaks)):
# 		#Bias? should this be randomized instead?
# 		tempScore = minScore(peaks[i], meanWordsList)
# 		if tempScore < minVal:
# 			minVal = tempScore
# 			seedIndex = i
# 	return (peaks.pop(seedIndex),peaks)

# The matrix array of characters is a list of [probA;probT;probG;probC] lists
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

#abstracts a peak seed into a mean
#CAREFUL, currently the only thing with action adding mean length
def abstract (peak):
	seq = peak[0]
	matrix = []
	for character in seq:
		matrix += [initProb(character)]
	return matrix

# picks a subsample from peaks, 
# then picks seed means from that by k++ w/o replacement
# ALSO, stop with no peaks in outliers - ?
def pickMeans (peaks, numMeans):
	#Now it deals with all len 0 clusters here
	if len(peaks) == 0:
		return []
	else:
		#New object with list() - not an alias
		nonReplacement = list(peaks)
		subSample = []
		subSampleSize = calcSubSize(numMeans, len(peaks))
		for x in xrange(subSampleSize):
			(thisPeak,rest) = sample(nonReplacement)
			nonReplacement = rest
			subSample += [thisPeak]
		#Now we create the list of new means
		seeds = []
		(peak, subSample) = sample(subSample)
		seed = abstract(peak)
		#Instantiate peakDistances so that the first alignment will always be better,
		#and replace the original value as the closest peak, so that the peak farthest from any is chosen next
		peakDistances = [0]*len(subSample)
		seeds += [seed]
		#PLEASE FIX numMeans guesses and index errors!!!!
		if numMeans > 1:
			for i in range(numMeans - 1):
				(seed, subSample, peakDistances) = kPlusPlus(seeds,subSample,peakDistances)
				seeds += [abstract(seed)]
		return seeds

#Assumes that the initial mean guess won't be larger that the number of peaks
#May implement so that it is the only function to use subsampling,
# assumes the initial k guess will reduce the number of sequences in any cluster significantly
def pickInitMeans (peaks, numMeans):\
	return pickMeans(peaks, numMeans)

# def pickInitMeans (peaks, numMeans):
# 	#New object with list() - not an alias
# 	nonReplacement = list(peaks)
# 	subSample = []
# 	subSampleSize = calcSubSize(numMeans, len(peaks))
# 	for x in xrange(subSampleSize):
# 		(thisPeak,rest) = sample(nonReplacement)
# 		nonReplacement = rest
# 		subSample += [thisPeak]
# 	return pickMeans(subSample, numMeans)

#Still assumes mean guess <= number of total peaks
def pickNewMeans (clusters, numMeans, clusterVariances):
	means = []
	outliers = clusters[0]
	numOutliers = len(outliers)
	#Perhaps divide the new seeds between outliers and old peaks based on size of outliers and variances?
	if numMeans > numOutliers:
		means += pickMeans(outliers, numOutliers)
		numMeans -= numOutliers
		#finds the highest variance cluster and it's index
		#!!!!!!!!!!!!!!!! Reevaluate the usefulness of this as opposed to a holistic sampling
		prevClusterVariances = list(enumerate(list(clusterVariances)))
		prevClusterVariances.sort(key = lambda seg: seg[1])
		while numMeans > 0:
			(j,maxVar) = prevClusterVariances.pop()
			cluster = clusters[j][1:]
			numPeaks = len(cluster)
			#variace of an empty mean is currently 0, but that could change,
			#so account for pulling from clusters without any means? - currently pickMeans returns [] in this case
			if numMeans > numPeaks:
				#If you pull apart an old mean, probs should throw out the old mean for the new ones.
				means += pickMeans(cluster, numPeaks)
				numMeans -= numPeaks
			else:
				means += pickMeans(clusters[j][1:], numMeans)
				numMeans -= numMeans
	else:
		means += pickMeans(outliers, numMeans)
	return means

