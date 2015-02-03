
# write individual peaks (or currently means)
def peakWriter(peak, outfile):
	sequence = peak[0]
	#line = footprint + '\t' + str(position) + '\t' + str(motif_start) + '\t' + str(seed) +'\n' (from input)
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
		clusterWriter(clusters[i], outfile)
	outfile.close()