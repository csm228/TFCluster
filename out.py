
def consensusCharacter(probArray):
	i = probArray.index(max(probArray))
	characters = ['A','T','G','C']
	return characters[i]

def consensusString(mean):
	consensus = ''
	for probArray in mean:
		consensus += consensusCharacter(probArray)
	return consensus

# write individual means
def meanWriter((mean,targetLength,targetIndex), specMotif, portion, outfile):
	line = consensusString(mean) + '\t' 'of mean alignment length ' + str(targetLength) + ' and mean alignment index ' + str(targetIndex) + ' with ' + specMotif + ' composing ' + str(portion) + '%\n'
	outfile.write(line)
	line = str(mean) + '\n'
	outfile.write(line)

# returns a sequence with only the motif uppercase
def withMotifUpper(seq,i1,motif):
	i2 = i1 + len(motif)
	newSeq = (seq[:i1]).lower() + seq[i1:i2] + (seq[i2:]).lower()
	return newSeq

# returns a sequence with only the motif lowercase
def withMotifLower(seq,i1,motif):
	i2 = i1 + len(motif)
	newSeq = seq[:i1] + (seq[i1:i2]).lower() + seq[i2:]
	return newSeq

# write individual peaks
def peakWriter(peak, outfile):
	sequence = peak[0]
	index = peak[2]
	motif = peak[3][1]
	newSeq = withMotifUpper(sequence,index,motif)
	line = newSeq + '\t' + motif + '\n'
	outfile.write(line)

# write clusters to output file
def clusterWriter(cluster, outfile):
	for peak in cluster:
		peakWriter(peak, outfile)
	outfile.write('\n')

def specificity(cluster):
	motifs = []
	for peak in cluster:
		#MEMOIZE maybe?
		motifs += [peak[3][1]]
	specMotif = 'no motifs'
	portion = 0
	numMotifs = len(motifs)
	if numMotifs > 0:
		specMotif = max(set(motifs), key=motifs.count)
		numOccurrences = motifs.count(specMotif)
		portion = (numOccurrences / float(numMotifs)) * 100.0
	return (specMotif,portion)


def writeOutput(filename, clusters):
	outfile = open(filename, 'w')
	outfile.write('Putative Outliers\n')
	clusterWriter(clusters[0], outfile)
	for i in range(1,len(clusters)):
		outfile.write('Cluster ' + str(i) + '\n')
		(specMotif, portion) = specificity(clusters[i][1:])
		meanWriter(clusters[i][0], specMotif, portion, outfile)
		clusterWriter(clusters[i][1:], outfile)
	outfile.close()