
import math

#The p-value to check against to stop the clustering algorithm. Change this
def probabilityThreshold = .05

#Use an asymptotic function with the Welch's t-test p-value as a seed for
#how far along the function to pick (closer to stable -> fewer new means)
#def calcSubSize (k,n):

def calcSubSize (numMeans, numPeaks):
	subSize = numPeaks/10 + 1
	return subSize
	
#Takes a random peak and removes it from the list of peaks,
#returns (peak, remaining peaks) as a tuple
def sample (peaks):
	index = random.randrange(len(peaks))
	return (peaks.pop(index),peaks)

def minScore (seq, meanWordsList):
	#align gives a tuple of index, score, so take the score
	minVal = score_of_align(seq, meanWordsList[0]))
	for meanWords in meanWordsList[1:]:
		minVal = min(minVal, score_of_align(seq,meanWords))
	return minVal

def kPlusPlus (means, peaks):
	meanWordsList = []
	for mean in means:
		meanWordsList += align.wordify(mean)
	seedIndex = 0
	minVal = minScore(peaks[0], meanWordsList)
	for i in range(1,len(peaks) - 1):
		#Bias? should this be randomized instead?
		if minScore(peaks[i], meanWordsList) < minVal:
			seedIndex = i
	return (peaks.pop(seedIndex),peaks)


# picks a subsample from peaks, 
# then picks seed means from that by k++ w/o replacement
# Need to implement what happens when numMeans > len(peaks)
def pickMeans (peaks, numMeans):
	#New object with list() - not an alias
	nonReplacement = list(peaks)
	subSample = []
	subSampleSize = calcSubSize(numMeans, len(peaks))
	for x in xrange(subSampleSize):
		(thisPeak,rest) = sample(nonReplacement)
		nonReplacement = rest
		subSample += thisPeak
	#Now we create the list of new means
	seeds = []
	(seed, subSample) = sample(subSample)
	seeds += [seed]
	for i in range(numMeans):
		()
		seeds += 
	return seeds


# try storing prev. cluster variance? - in main
def variance (cluster):
	sumVar = 0
	for seq in cluster[1:]
		(index, score) = align(seq, cluster[0])
		sumVar += score
	#Integer division?
	return (sumVar / (len(cluster) - 1))

# returns the standard deviation of a cluster
# measured by alignment score from the mean
def stdDev (cluster):
	return math.sqrt(variance(cluster))

def parse(file):

#Tie this to subsample
def guessInitMeans(peaks):
	return 5

def welchTest (prevClusters,currClusters):
	prevMeans = []
	currMeans = []
	for cluster in prevClusters[1:]:
		prevMeans += [cluster[0]]
	for cluster in currClusters[1:]:
		currMeans += [cluster[0]]
	(t_stat,p_val) = scipy.stats.ttest_ind(prevMeans, currMeans, equal_var = False)
	return p_val

# print a line when the variance goes back up between means
# Requires a list of peaks where the sequence is the first element of the peak
def main (filename):
	peaks = parse(filename)
	prevClusters = []
	#The extra list on the end is for outliers
	currClusters = [[]] + pickMeans(peaks,guessInitMeans(peaks))
	while welchTest(prevClusters,currClusters) > probabilityThreshold:
		cluster(peaks,means)
