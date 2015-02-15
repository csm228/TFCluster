import align

#Currently we just throw out peaks that the alignment algorithm
#doesn't align (returns 0 instead)
outlierThreshold = 0

allocationCessationThreshold = 1000

#In case the datatype changes, or we want to
#remember the cluster of the peak by storing it in the peak
def group (peak, meanIndex, clusters):
	clusters[meanIndex] += [peak]

#Generates blank prototypes. Shouldn't be affected by paring
def initializePrototypes(means):
	prototypes = []
	for m in range(len(means)):
		#The prototypical mean, generated during allocation
		prototype = []
		for n in range(len(means[m])):
			prototype += [[0,0,0,0]]
		prototypes += [prototype]
	return prototypes

#Counts the character's contribution to the mean 
# at that position in the sequence and mean
def count (character, probArray):
	if character == 'A':
		probArray[0] += 1
	elif character == 'T':
		probArray[1] += 1
	elif character == 'G':
		probArray[2] += 1
	elif character == 'C':
		probArray[3] += 1
	else:
		print "Incorrect string in peaks (recentering), implement error handling"
	return

def account(peak, mean, i):
	for j in range(max(0,-i),min(len(peak[0]), meanLength - i)):
		print prototype[j+i]
		count(peak[0][j],prototype[j+i])
		print prototype[j+i]

#The allocation step of k-means clustering
#Now also adds to the prototypical means for the recentering step
def allocate (peaks, means, alignmentMatrix):
	#Outliers are stored in the first cluster (list) in clusters
	clusters = [[]]
	prototypes = initializePrototypes(means)
	for i in range(len(peaks)):
		#the score matrix needs to start with a value for the outlier cluster
		alignments = [(0,0)] + alignmentMatrix[i]
		for j in range(len(means)):
			(i,score) = alignmentMatrix[i][j]
			scores += [score]
		maxScore = scores[0]
		#the closest mean, instantiated as the outlier cluster
		nearest = 0
		for k in range(1,len(scores)):
			#How to resolve ties?
			if scores[k] > maxScore:
				maxScore = scores[k]
				nearest = k
		if maxScore <= outlierThreshold:
			group(peaks[i], 0, clusters)
		else:
			group(peaks[i], nearest, clusters)
	return clusters

# def allocate (peaks, means, alignmentMatrix):
# 	#Outliers are stored in the first cluster (list) in clusters
# 	clusters = [[]]
# 	meanWordLists = []
# 	for mean in means:
# 		clusters += [[mean]]
# 		meanWordLists += [align.wordify(mean)]
# 	for peak in peaks:
# 		#the closest mean, instantiated as the outlier cluster
# 		scores = []
# 		for meanWords in meanWordLists:
# 			# (i,score) = align.align(peak,meanWords)
# 			scores += [align.score_of_align(peak,meanWords)]
# 		maxScore = scores[0]
# 		nearest = 0
# 		for j in range(1,len(scores)):
# 			#How to resolve ties?
# 			if scores[j] > maxScore:
# 				maxScore = scores[j]
# 				nearest = j
# 		if maxScore <= outlierThreshold:
# 			group(peak, 0, clusters)
# 		else:
# 			group(peak, nearest, clusters)
# 	return clusters


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
def recenter (clusters, deltaMeans, alignmentMatrix):
	for j in range(1,len(clusters)):
		#The prototypical mean, generated at recentering
		prototype = []
		cluster = clusters[j]
		#cluster[0] is the previous mean, j = 0 is the outlier cluster
		for i in range(len(cluster[0])):
			prototype += [[0,0,0,0]]
		#calculating the distribution in bases of the mean
		meanWords = align.wordify(cluster[0])
		meanLength = len(cluster[0])
		for peak in cluster[1:]:
			#MEMOIZE THE HECK OUT OF THIS PLEASE
			i = align.index_of_align(peak,meanWords)
			#May want to change this for extending means - when peaks flow over
			#Currently, it just keeps the original seed length & alignment
			#ALSO, accounts for alignments prior to the beginning of the mean,
			# or flowing over the end
			for j in range(max(0,-i),min(len(peak[0]), meanLength - i)):
				print prototype[j+i]
				count(peak[0][j],prototype[j+i])
				print prototype[j+i]
		for loc in prototype:
			total = 0
			for prob in loc:
				# print prob
				total += prob
			#should there ever be a mean without peaks? I don't think so
			if total != 0:
				for prob in loc:
					# print total
					prob /= total
		#Here is where highly variant means should be thrown out,
		#but need to allow for the first run with a mean - 
		# the seed will always have high change
		deltaMeans += difference(cluster[0],prototype)
		cluster[0] = prototype
	return deltaMeans

# def recenter (clusters,deltaMeans):
# 	for cluster in clusters[1:]:
# 		#The prototypical mean, generated at 
# 		prototype = []
# 		#cluster[0] is the previous mean
# 		for i in range(len(cluster[0])):
# 			prototype += [[0,0,0,0]]
# 		#calculating the distribution in bases of the mean
# 		meanWords = align.wordify(cluster[0])
# 		meanLength = len(cluster[0])
# 		for peak in cluster[1:]:
# 			#MEMOIZE THE HECK OUT OF THIS PLEASE
# 			i = align.index_of_align(peak,meanWords)
# 			#May want to change this for extending means - when peaks flow over
# 			#Currently, it just keeps the original seed length & alignment
# 			#ALSO, accounts for alignments prior to the beginning of the mean,
# 			# or flowing over the end
# 			for j in range(max(0,-i),min(len(peak[0]), meanLength - i)):
# 				print prototype[j+i]
# 				count(peak[0][j],prototype[j+i])
# 				print prototype[j+i]
# 		for loc in prototype:
# 			total = 0
# 			for prob in loc:
# 				# print prob
# 				total += prob
# 			#should there ever be a mean without peaks? I don't think so
# 			if total != 0:
# 				for prob in loc:
# 					# print total
# 					prob /= total
# 		#Here is where highly variant means should be thrown out,
# 		#but need to allow for the first run with a mean - 
# 		# the seed will always have high change
# 		deltaMeans += difference(cluster[0],prototype)
# 		cluster[0] = prototype
# 	return deltaMeans

#to extract the means so that they can be given back to main
def extractMeans (clusters):
	means = []
	for cluster in clusters[1:]:
		means += [cluster[0]]
	return means

#The k-means clustering algorithm, managing the termination of clustering under a guessed number of means
def cluster (peaks, means):
	#The total change in the means by alignment score, this should be eventually replaced
	deltaMeans = 1001
	initDeltaMeans = 0
	clusters = []
	# print means[0]
	while deltaMeans > allocationCessationThreshold:
		alignmentMatrix = generate_align_matrix(peaks,means)
		clusters = allocate(peaks,means,alignmentMatrix)
		# print clusters[1][0]
		deltaMeans = recenter(clusters,initDeltaMeans,alignmentMatrix)
		means = extractMeans(clusters)
	return (extractMeans(clusters),clusters)