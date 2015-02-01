
#how long the words will be
def wordLength = 5

#Lower bound for inclusion in high scoring words
def highScoreThreshold = 3.0

#Distance that two words must be under to group into a segment pair
def segPairWordMaxDist = 8

#The matrix array of characters is 
def compare (character, probArray):
	if character == A:
		return probArray[0]
	if character == T:
		return probArray[1]
	if character == G:
		return probArray[2]
	if character == C:
		return probArray[3]
	else print "Incorrect string in peaks, implement error handling"

def scorePair (seqWord,meanWord):
	#if misalignments become an issue, use weighted mean?
	score = 0
	for i in range(wordLength):
		score += compare(seqWord[i], meanWord[i])
	return score

#Takes a sequence (or mean) and returns the list of words of length wordLength
#Used in the alignment function
def wordify (seq):
	seqWords = []
	for i in range((len(seq)-wordLength)):
		seqWords += seq[i:(i+wordLength)]
	return seqWords

#Generates an index of alignment and an alignment score
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
	if len(highScoreWords) > 0:
		segmentPairs = []
		for k in range(len(highScoreWords)):
			(score1,i1,j1) = highScoreWords[k]
			for l in range(len(highScoreWords) - k):
				(score2,i2,j2) = highScoreWords[k]
				#Second high-scoring word must be close enough, on the same alignment
				dist = abs(i1-i2)
				if (dist < segPairWordMaxDist) && (i1 - i2 + j2 = j1):
					segScore = 0
					for m in range(dist):
						seqScore += scoreMatrix[min(i1,i2)+m][min(j1,j2)+m]
					segmentPairs += [(min(j1,j2)-min(i1,i2),seqScore)]
		#FIX THIS - sort function? implement search trees for ^ ?
		if len(segmentPairs) > 0:
			segmentPairs.sort(key = lambda seg: seg[1], reverse=True)
			return segmentPairs[0]
		highScoreWords.sort(reverse=True)
		(score,i,j) = highScoreWords[0]
		return (j-i,score)
	return (0,0)


def align_mean (peak, mean):
	return align(peak, wordify(mean))

#kinda sad functions, just left them in, used in cluster.py
def index_of_align(seq, meanWords):
	(i,score) = align(seq,meanWords)
	return i

def score_of_align(seq, meanWords):
	(i,score) = align(seq,meanWords)
	return score