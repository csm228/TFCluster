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

def kPlusPlus (means, peaks):
	print "kPlusPlus"
	meanWordsList = []
	for mean in means:
		meanWordsList += [align.wordify(mean)]
	seedIndex = 0
	minVal = minScore(peaks[0], meanWordsList)
	for i in range(1,len(peaks)):
		#Bias? should this be randomized instead?
		tempScore = minScore(peaks[i], meanWordsList)
		if tempScore < minVal:
			minVal = tempScore
			seedIndex = i
	return (peaks.pop(seedIndex),peaks)

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

# def initProb (character):
# 	if character == 'A':
# 		return [1,0,0,0]
# 	elif character == 'T':
# 		return [0,1,0,0]
# 	elif character == 'G':
# 		return [0,0,1,0]
# 	elif character == 'C':
# 		return [0,0,0,1]
# 	else:
# 		print "Incorrect string in peaks, implement error handling"

#This rearrangement of seed generation may reduce mean locations with [0,0,0,0]
#OR NOT :[
# def initProb (character):
# 	if character == 'A':
# 		return [5,1,1,1]
# 	elif character == 'T':
# 		return [1,5,1,1]
# 	elif character == 'G':
# 		return [1,1,5,1]
# 	elif character == 'C':
# 		return [1,1,1,5]
# 	print "Incorrect string in peaks, implement error handling"

#abstracts a peak seed into a mean
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
		(seed, subSample) = sample(subSample)
		seeds += [abstract(seed)]
		#PLEASE FIX numMeans guesses and index errors!!!!
		if numMeans > 1:
			for i in range(numMeans - 1):
				(seed, subSample) = kPlusPlus(seeds,subSample)
				seeds += [abstract(seed)]
		return seeds

#Assumes that the initial mean guess won't be larger that the number of peaks
#May implement so that it is the only function to use subsampling,
# assumes the initial k guess will reduce the number of sequences in any cluster significantly
def pickInitMeans (peaks, numMeans):
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

