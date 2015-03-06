
#how long the words will be
wordLength = 5

#Lower bound for inclusion in high scoring words
highScoreThreshold = 4

#Distance that two words must be under to group into a segment pair
segPairWordMaxDist = 11

#The matrix array of characters is a list of [probA;probT;probG;probC] lists
def compare (character, probArray):
	if character == 'A':
		return probArray[0]
	elif character == 'T':
		return probArray[1]
	elif character == 'G':
		return probArray[2]
	elif character == 'C':
		return probArray[3]
	else:
		print "Incorrect string in peaks, implement error handling"

#assume both have the same length, generally = wordLength
def scorePair (seqWord,meanWord):
	#if misalignments become an issue, use weighted mean?
	score = 0
	for i in range(len(seqWord)):
		score += compare(seqWord[i], meanWord[i])
	return score

#Takes a sequence (or mean) and returns the list of words of length wordLength
#Used in the alignment function
def wordify (seq):
	seqWords = []
	for i in range((len(seq)-wordLength)):
		seqWords += [seq[i:(i+wordLength)]]
	return seqWords

#Generates an index of alignment and an alignment score
# in the form (index, score)
# where the index is the displacement of the peak sequence from the mean
def align (peak, meanWords):
	seqWords = wordify(peak[0])
	scoreMatrix = []
	highScoreWords = []
	for i in range(len(seqWords)):
		scores = []
		for j in range(len(meanWords)):
			score = scorePair(seqWords[i],meanWords[j])
			scores += [score]
			if score > highScoreThreshold:
				highScoreWords += [(score,i,j)]
		scoreMatrix += [scores]
	# print highScoreWords[0]
	# print scoreMatrix
	if len(highScoreWords) > 0:
		segmentPairs = []
		for k in range(len(highScoreWords)):
			(score1,i1,j1) = highScoreWords[k]
			for l in range(k + 1,len(highScoreWords)):
				# print highScoreWords[k]
				# print highScoreWords[l]
				# print '\n'
				(score2,i2,j2) = highScoreWords[l]
				#Second high-scoring word must be close enough, on the same alignment
				dist = abs(i1-i2)
				if (dist < segPairWordMaxDist) and (i1 - i2 + j2 == j1):
					segScore = 0
					# for m in range(dist):
					# 	segScore += scoreMatrix[min(i1,i2)+m][min(j1,j2)+m]
					peakSeg = peak[0][i1:i2+wordLength]
					meanSeg = []
					#This could reaaaally use MEMOIZATION
					for j in range(j1,j2):
						meanSeg += [meanWords[j][0]]
					meanSeg += meanWords[j2]
					#So now, doesn't double-count characters in the center
					segScore = scorePair(peakSeg,meanSeg)
					segmentPairs += [(min(j1,j2)-min(i1,i2),segScore)]
		#FIX THIS - sort function? implement search trees for ^ ?
		if len(segmentPairs) > 0:
			segmentPairs.sort(key = lambda seg: seg[1], reverse=True)
			return segmentPairs[0]
		highScoreWords.sort(reverse=True)
		(score,i,j) = highScoreWords[0]
		return (j-i,score)
	return (0,0)


# def align_mean (peak, mean):
# 	return align(peak, wordify(mean))

#kinda sad functions, index is only used in recentering,
# and only the score is used everywhere else: so separate them entirely?
def index_of_align(peak, meanWords):
	(i,score) = align(peak,meanWords)
	return i

def score_of_align(peak, meanWords):
	(i,score) = align(peak,meanWords)
	return score

def generate_align_matrix(peaks,means):
	meanWordsList = []
	for mean in means:
		meanWordsList += [wordify(mean)]
	alignmentMatrix = []
	for i in range(len(peaks)):
		peak = peaks[i]
		peakAlignments = []
		for j in range(len(means)):
			peakAlignments += [align(peak,meanWordsList[j])]
		alignmentMatrix += [peakAlignments]
	return alignmentMatrix