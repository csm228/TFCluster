# make footprints, write footprints to output file
def fprintWriter(matrix, length, type, outfile, position, seed):
	footprint = ""	
	added_bases = length - len(matrix)
	motif_start = random.randint(0, added_bases-1)
	motif_end = added_bases - motif_start
	for i in range(motif_start):
		base = sampleDistr(RANDOM_SEQ)
		footprint = footprint + base		
	for distr in matrix:
		base = sampleDistr(distr)
		footprint = footprint + base		
	for i in range(motif_end):
		base = sampleDistr(RANDOM_SEQ)
		footprint = footprint + base
	line = footprint + '\t' + str(position) + '\t' + str(motif_start) + '\t' + str(seed) +'\n'
	outfile.write(line)

def writeOutput(filename, clusters)
outfile = open(sys.argv[5], 'w')
for i in range(f_count):
	use_seed = random.choice(seeds)
	foot_length = random.randint(foot_range[0],foot_range[1])
	fprintWriter(use_seed[2], foot_length, use_seed[0], outfile, i, use_seed)
outfile.close()