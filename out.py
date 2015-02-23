
def consensusCharacter(probArray):
	i = probArray.index(max(probArray))
	characters = ['A','T','G','C']
	return characters[i]

def consensusString(mean):
	consensus = ''
	for probArray in mean:
		consensus += consensusCharacter(probArray)
	return consensus

# write individual peaks (or currently means)
def meanWriter(mean, outfile):
	line = consensusString(mean) + '\n'
	outfile.write(line)
	line = str(mean) + '\n'
	outfile.write(line)

# returns a sequence with only the motif uppercase
def withMotif(seq,i1,motif):
	i2 = i1 + len(motif)
	newSeq = (seq[:i1]).lower() + ' ' + seq[i1:i2] + ' ' + (seq[i2:]).lower()
	return newSeq

# write individual peaks
def peakWriter(peak, outfile):
	sequence = peak[0]
	index = peak[2]
	motif = peak[3][1]
	newSeq = withMotif(sequence,index,motif)
	line = newSeq + '\t' + motif + '\n'
	outfile.write(line)

# write clusters to output file
def clusterWriter(cluster, outfile):
	for peak in cluster:
		peakWriter(peak, outfile)
	outfile.write('\n')

def writeOutput(filename, clusters):
	outfile = open(filename, 'w')
	outfile.write('Putative Outliers\n')
	clusterWriter(clusters[0], outfile)
	for i in range(1,len(clusters)):
		outfile.write('Cluster ' + str(i) + '\n')
		meanWriter(clusters[i][0], outfile)
		clusterWriter(clusters[i][1:], outfile)
	outfile.close()