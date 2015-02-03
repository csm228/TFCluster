import align

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
		meanWordLists += [align.wordify(mean)]
	for peak in peaks:
		#the closest mean, instantiated as the outlier cluster
		scores = []
		for mean in means:
			(i,score) = align.align(peak,mean)
			scores += score
		maxScore = scores[0]
		nearest = 0
		for j in range(1,len(scores)):
			#How to resolve ties?
			if scores[j] > maxScore:
				maxScore = scores[j]
				nearest = j
		if maxScore <= outlierThreshold:
			group(peak, 0, clusters)
		else:
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
	else:
		print "Incorrect string in peaks (recentering), implement error handling"
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
def recenter (clusters,deltaMeans):
	for cluster in clusters[1:]:
		#The prototypical mean, generated at 
		prototype = []
		#cluster[0] is the previous mean
		for i in range(len(cluster[0])):
			prototype += [[0,0,0,0]]
		#calculating the distribution in bases of the mean
		meanWords = align.wordify(cluster[0])
		for peak in cluster:
			i = align.index_of_align(peak,meanWords)
			#May want to change this for extending means - when peaks flow over
			#Currently, it just keeps the original seed length & alignment
			for j in range(min(len(peak),(len(cluster[0]) - i))):
				count(peak[0][j],prototype[i+j])
		for loc in prototype:
			count = 0
			for prob in loc:
				count += prob
			for prob in loc:
				prob /= count
		#Here is where highly variant means should be thrown out,
		#but need to allow for the first run with a mean - 
		# the seed will always have high change
		deltaMeans += difference(cluster[0],prototype)
		cluster[0] = prototype
		return deltaMeans

#to extract the means so that they can be given back to main
def extractMeans (clusters):
	means = []
	for cluster in clusters[1:]:
		means += cluster[0]
	return means

#The k-means clustering algorithm, managing the termination of the 
def cluster (peaks, means):
	#The total change in the means by alignment score, this should be eventually replaced
	deltaMeans = 1001
	clusters = []
	while deltaMeans > allocationCessationThreshold:
		clusters = allocate (peaks,means)
		deltaMeans = recenter (clusters,0)
	return (extractMeans(clusters),clusters)