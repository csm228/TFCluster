import align
import clustering
import seed
import paring

import random
import math
import copy #need to make deep copies anywhere?

import scipy
from scipy import stats


#The p-value to check against to stop the clustering algorithm. Change this
probabilityThreshold = .6

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
		print "Incorrect string in peaks"

#abstracts a peak seed into a mean
#CAREFUL, currently the only thing with action adding mean length
def abstract (peak):
	seq = peak[0]
	matrix = []
	for character in seq:
		matrix += [initProb(character)]
	return (matrix, -1)


def pickTwo(peaks):
	numPeaks = len(peaks)
	if numPeaks > 1:
		index = random.randrange(numPeaks)
		seed1 = abstract(peaks[index])
		meanWords = align.wordify(seed1)
		#Pick randomly? Take one that aligns something else? (one-off)
		farthest = 0
		farthestScore = align.score_of_align(peaks[0],meanWords,-1)
		for i in range(1,numPeaks):
			#score_of_align needs targetLength now, -1 for inconsequential
			score = align.score_of_align(peaks[i],meanWords,-1)
			if score < farthestScore:
				farthestScore = score
				farthest = i
		seed2 = abstract(peaks[farthest])
		return [seed1, seed2]
	# elif numPeaks == 1:
	# 	return [abstract(peaks[0])]
	else:
		return []



#Now using varAlignmentMatrix
def variance (cluster,meanNum,varAlignmentMatrix):
	sumVar = 0
	meanLength = len(cluster[0])
	for (i,m,score,length) in varAlignmentMatrix[meanNum]:
		#This definition of distance is kinda inconsistent -
		# alignment gives word and segment scores, not whole sequence alignments
		sumVar += (meanLength - score)**2
	n = len(cluster) - 1
	if n != len(varAlignmentMatrix[meanNum]):
		print "the number of peaks in the cluster doesn\'t match the number of alignments - invariant invalidity"
	if n > 0:
		#Integer division is "//" while "/" does floating point
		return (sumVar / n)	
	#What do you do if NOTHING gets clustered with the mean?
	else:
		return 0

# def variance (cluster,meanNum,alignmentMatrix):
# 	sumVar = 0
# 	meanLength = len(cluster[0])
# 	for peak in cluster[1:]:
# 		peakNum = peak[1]
# 		(i,m,score,length) = alignmentMatrix[peakNum][meanNum]
# 		#This definition of distance is kinda inconsistent -
# 		# alignment gives word and segment scores, not whole sequence alignments
# 		sumVar += (meanLength - score)**2
# 	n = len(cluster) - 1
# 	if n > 0:
# 		#Integer division is "//" while "/" does floating point
# 		return (sumVar / n)	
# 	#What do you do if NOTHING gets clustered with the mean?
# 	else:
# 		return 0

# returns the standard deviation of a cluster
# measured by alignment score from the mean
def stdDev (cluster):
	return math.sqrt(variance(cluster))

#New with lower momoization costs (only visible in variance calculation)
def welchTest (currClusters,alignmentMatrix,prevClusterVariances):
	currClusterVariances = [0]
	for c in range(len(currClusters)-1):
		clusterVariance = variance(currClusters[c+1],c,alignmentMatrix)
		if clusterVariance > 0:
			currClusterVariances += [clusterVariance]
	# # What if nothing is in either cluster, all in outliers?
	# # What if this happens perpetually?
	# if currClusterVariances == []:
	# 	return (0,[0])
	(t_stat,p_val) = scipy.stats.ttest_ind(prevClusterVariances, currClusterVariances, equal_var = False)
	print p_val
	return (p_val,currClusterVariances)

def alignmentWelchTest(currClusters,varAlignmentMatrix,prevClusterScores):
	currClusterScores1 = []
	print varAlignmentMatrix[0]
	for (i,m,score,length) in varAlignmentMatrix[0]:
		currClusterScores1 += [score]
	print currClusterScores1
	(t_stat1,p_val1) = scipy.stats.ttest_ind(prevClusterScores, currClusterScores1, equal_var = False)
	print varAlignmentMatrix[1]
	currClusterScores2 = []
	for (i,m,score,length) in varAlignmentMatrix[1]:
		currClusterScores2 += [score]
	print currClusterScores2
	(t_stat2,p_val2) = scipy.stats.ttest_ind(prevClusterScores, currClusterScores2, equal_var = False)
	p_val = (p_val1 * p_val2)
	print p_val
	return (p_val,[currClusterScores1,currClusterScores2])

#temporary initial step for the welch's test (?)
def clustrifyMeans (means):
	listifiedMeans = []
	for i in range(len(means)):
		#Double brackets, as single brackets just recostructs means
		listifiedMeans += [[means[i]]]
	# print listifiedMeans[0]
	return listifiedMeans

def bifurcate (cluster, prevClusterScore):
	numPeaks = len(cluster) - 1
	if numPeaks == 1:
		return [[abstract(peaks[0]),peaks[0]]]
	if numPeaks == 0:
		return []
	peaks = cluster[1:]
	means = pickTwo(peaks)
	#The extra list at the beginning is for outliers,and is initialized with all peaks
	clusters = [peaks] + clustrifyMeans(means)
	alignmentMatrix = align.generate_align_matrix(peaks,means)
	(means,clusters) = clustering.cluster(peaks,means,alignmentMatrix)
	varAlignmentMatrix = align.generate_var_align_matrix(clusters)
	print 'finished a binary branching'
	# (p_val, clusterVariances) = welchTest(clusters,varAlignmentMatrix,[prevClusterVariance])
	(p_val, clusterScores) = alignmentWelchTest(clusters,varAlignmentMatrix,[prevClusterScore])
	if p_val > probabilityThreshold:
		return [cluster]
	else:
		binaryClusters = main(clusters[0]) + bifurcate(clusters[1],clusterScores[0]) + bifurcate(clusters[2],clusterScores[1])
		return binaryClusters

#SO MUCH MEMOIZATION
# print a line when the variance goes back up between means
# Requires a list of peaks where the sequence is the first element of the peak
def main (peaks):
	numPeaks = len(peaks)
	if numPeaks == 1:
		return [[abstract(peaks[0]),peaks[0]]]
	if numPeaks == 0:
		return []
	means = pickTwo(peaks)
	print peaks[0]
	#The extra list at the beginning is for outliers,and is initialized with all peaks
	clusters = [peaks] + clustrifyMeans(means)
	alignmentMatrix = align.generate_align_matrix(peaks,means)
	prevClusterScore = [0]*numPeaks #just something so that the first Welch's test doesn't cause termination
	(means,clusters) = clustering.cluster(peaks,means,alignmentMatrix)
	varAlignmentMatrix = align.generate_var_align_matrix(clusters)
	print 'finished an outlier clustering (main)'
	# (p_val, clusterVariances) = welchTest(clusters,varAlignmentMatrix,[prevClusterVariance])
	(p_val, clusterScores) = alignmentWelchTest(clusters,varAlignmentMatrix,[prevClusterScore])
	# if p_val > probabilityThreshold:
	# 	means = paring.paredMeans(means,varAlignmentMatrix)
	# 	numNewMeans = guessNewMeans(peaks, means, p_val)
	# 	#currently, no correlation between how many means duplicated/dropped in paring
	# 	#and how many and from where they are added in mean picking
	# 	means += seed.pickNewMeans(clusters, numNewMeans, clusterVariances)
	# 	alignmentMatrix = align.generate_align_matrix(peaks,means)
	# 	(means,clusters) = cluster.cluster(peaks,means,alignmentMatrix)
	# 	varAlignmentMatrix = align.generate_var_align_matrix(clusters)
	# 	(p_val, clusterVariances) = welchTest(clusters,varAlignmentMatrix,clusterVariances)
	# 	print 'finished clustering of subsequent k guess'
	# 	n += 1
	binaryClusters = bifurcate(clusters[0],clusterScores[0]) + bifurcate(clusters[1],clusterScores[1])
	return binaryClusters