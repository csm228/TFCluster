import align
import cluster
import seed
import paring

import random
import math
import copy #need to make deep copies anywhere?

import scipy
from scipy import stats


#The p-value to check against to stop the clustering algorithm. Change this
probabilityThreshold = .96


# try storing prev. cluster variance? - in main
def variance (cluster,j,alignmentMatrix):
	sumVar = 0
	meanLength = len(cluster[0])
	for peak in cluster[1:]:
		#This definition of distance is kinda inconsistent -
		# alignment gives word and segment scores, not whole sequence alignments
		sumVar += (meanLength - alignmentMatrix[peak[1]][j])**2
	n = len(cluster) - 1
	if n > 0:
		#Integer division is "//" while "/" does floating point
		return (sumVar / n)	
	#What do you do if NOTHING gets clustered with the mean?
	else:
		return 0

# returns the standard deviation of a cluster
# measured by alignment score from the mean
def stdDev (cluster):
	return math.sqrt(variance(cluster))

#Use an asymptotic function with the Welch's t-test p-value as a seed for
#how far along the function to pick (closer to stable -> fewer new means)

#Tie this to subsample in the future
def guessInitMeans(peaks):
	return (len(peaks) // 15) + 5

#Guesses how many more means will be needed
#Just iterates by 5. Ouch.
def guessNewMeans(peaks, means, p_val):
	return 5

def welchTest (prevClusters,currClusters,currAlignmentMatrix,prevAlignmentMatrix):
	prevClusterVariances = []
	currClusterVariances = []
	for c1 in range(1,len(prevClusters)):
		prevClusterVariances += [variance(prevClusters[c1],c1,prevAlignmentMatrix)]
	for c2 in range(1,len(currClusters)):
		currClusterVariances += [variance(currCluster[c2],c2,currAlignmentMatrix)]
	(t_stat,p_val) = scipy.stats.ttest_ind(prevClusterVariances, currClusterVariances, equal_var = False)
	print p_val
	return p_val

#temporary initial step for the welch's test (?)
def clustrifyMeans (means):
	listifiedMeans = []
	for i in range(len(means)):
		#Double brackets, as single brackets just recostructs means
		listifiedMeans += [[means[i]]]
	# print listifiedMeans[0]
	return listifiedMeans

#SO MUCH MEMOIZATION
# print a line when the variance goes back up between means
# Requires a list of peaks where the sequence is the first element of the peak
def main (peaks):
	prevClusters = []
	means = seed.pickInitMeans(peaks,guessInitMeans(peaks))
	print peaks[0]
	#The extra list at the beginning is for outliers,and is initialized with all peaks
	currClusters = [peaks] + clustrifyMeans(means)
	prevClusters = currClusters
	alignmentMatrix = align.generate_align_matrix(peaks,means)
	clusterVariances = [0]
	print 'first runthrough of clustering'
	(means,currClusters) = cluster.cluster(peaks,means,currAlignmentMatrix)
	currAlignmentMatrix = align.generate_align_matrix(peaks,means)
	p_val = welchTest(prevClusters,currClusters,currAlignmentMatrix,prevAlignmentMatrix)
	print 'starting welch\'s t-test clustering with centroid means'
	while p_val < probabilityThreshold:
		means = paring.paredMeans(means)
		numNewMeans = guessNewMeans(peaks, means, p_val)
		means += pickMeans(currClusters, numNewMeans)
		prevClusters = currClusters
		prevAlignmentMatrix = currAlignmentMatrix
		currAlignmentMatrix = align.generate_align_matrix(peaks,means)
		(means,currClusters) = cluster.cluster(peaks,means,currAlignmentMatrix)
		p_val = welchTest(prevClusters,currClusters,currAlignmentMatrix,prevAlignmentMatrix)
		print 'finished clustering of subsequent k guess'
	return currClusters