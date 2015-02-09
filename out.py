
# write individual peaks (or currently means)
def meanWriter(mean, outfile):
	line = str(mean) + '\n'
	outfile.write(line)

# write individual peaks
def peakWriter(peak, outfile):
	sequence = peak[0]
	line = sequence + '\n'
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