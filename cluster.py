
#Currently we just throw out peaks that the alignment algorithm
#doesn't align (returns 0 instead)
outlierThreshold = 0

allocationCessationThreshold = 1000

#In case the datatype changes, or we want to
#remember the cluster of the peak by storing it in the peak
def group (peak, meanIndex, clusters):
	clusters[meanIndex] += [peak]


#The allocation step of k-means clustering
def allocate (peaks, means):
	#Outliers are stored in the first cluster (list) in clusters
	clusters = [[]]
	meanWordLists = []
	for mean in means:
		clusters += [[mean]]
		meanWordLists += [Align.wordify(mean)]
	for peak in peaks:
		#the closest mean, instantiated as the outlier cluster
		scores = []
		for j in range(len(means)):
			(i,score) = Align.align(peak,mean)
			scores += 
		if max_score <= outlierThreshold:
			group(peak, len(clusters) - 1, clusters)
		group(peak, nearest, clusters)
	return clusters

#Counts the character's contribution to the mean 
# at that position in the sequence and mean
def count (character, probArray):
	if character == A:
		probArray[0] += 1
	if character == T:
		probArray[0] += 1
	if character == G:
		probArray[0] += 1
	if character == C:
		probArray[0] += 1
	else print "Incorrect string in peaks (recentering), implement error handling"
	return

#Reliant on means of the same size....
def difference (prevMean, currMean):
	diff = 0
	for i in range(len(currMean)):
		#This should only go to four - len(currMean[i])
		for j in range(4):
			diff += abs(prevMean[i][j] - currMean[i][j])
	return diff


#if recentered mean moves from seed (previous mean), throw out and pick only from outliers
#recenter on the fly: store number in cluster and just add then average??
#The calculation step of the "centroids" of the clusters
def recenter (clusters,deltaMean):
	for cluster in clusters[1:]:
		#The prototypical mean, generated at 
		prototype = []
		#cluster[0] is the previous mean
		for i in range(len(cluster[0])):
			prototype += [[0,0,0,0]]
		#calculating the distribution in bases of the mean
		meanWords = Align.wordify(cluster[0])
		for peak in cluster:
			i = index_of_align(peak,meanWords)
			#May want to change this for extending means - when peaks flow over
			#Currently, it just keeps the original seed length & alignment
			for j in range(min(len(peak),len(cluster[0] - i))):
				count(peak[0][j],prototype[i+j])
		for loc in prototype:
			count = 0
			for prob in loc:
				count += prob
			for prob in loc:
				prob /= count
		deltaMean += difference(cluster[0],prototype)
		cluster[0] = prototype
		return delta_mean

#The k-means clustering algorithm, managing the termination of the 
def cluster (peaks, means):
	#The total change in the means by alignment score
	deltaMean = 1001
    clusters = []
	while delta_mean > allocationCessationThreshold:
		clusters = allocate (peaks,means)
		deltaMean = recenter (peaks,means,clusters)
	return clusters