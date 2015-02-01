# 
# USE: python motifMaker.py <number of seeds> <number of footprints> <min footprint length> <max footprint length> <outfile name>
# OUTPUT: Footprint sequence, sequence position(number), motif start position, Seed info (list)
import random
import sys

RANDOM_SEQ = {'A': 25, 'T': 50, 'C': 75, 'G': 100}

# Convert distributions to cumulative distributions
def convertDistr(matrix):
	for distr in matrix:
		A, T, C, G = distr['A'], distr['T'], distr['C'], distr['G'] 
		if (A + T + C + G) > 100:
			continue
		else:	
			distr['A'] = A
			distr['T'] = A + T
			distr['C'] = A + T + C
			distr['G'] = A + T + C + G
	return matrix

# Make footprint motifs	
def makeMotif(distr, base_conserv, min_length, max_length, split):	
	motif = ""
	matrix = []
	motif_length = random.randint(min_length,max_length)	
	for i in range(motif_length):
		sample = random.randint(1,100)
		if sample <= distr.get('A'):
			motif += 'A'
		elif sample <= distr.get('T'):
			motif += 'T'
		elif sample <= distr.get('C'):
			motif += 'C'
		else:
			motif += 'G'	
	if split == 0:
		for base in motif:
			bases = ['A','T','C','G']
			base_matrix = {}
			base_prob = random.randint(base_conserv, 100)
			base_matrix.update({base:base_prob})
			bases.remove(base)
			for i in range(2):
				new_base = random.choice(bases)
				prob_left = 100 - base_prob
				new_prob = random.randint(0, prob_left)
				base_matrix.update({new_base:new_prob})
				base_prob = base_prob + new_prob
				bases.remove(new_base)
			new_base = random.choice(bases)
			prob_left = 100-base_prob
			base_matrix.update({new_base:prob_left})
			matrix.append(base_matrix)
	else:
		mid_start = random.randint(4,(motif_length-4))
		mid_length = random.randint(1,(motif_length/3))
		start_motif = motif[0:(mid_start-1)]
		mid_motif = motif[(mid_start-1):(mid_start + mid_length)]
		end_motif = motif[(mid_start + mid_length):motif_length]
		for base in start_motif:
			bases = ['A','T','C','G']
			base_matrix = {}
			base_prob = random.randint(base_conserv, 100)
			base_matrix.update({base:base_prob})
			bases.remove(base)
			for i in range(2):
				new_base = random.choice(bases)
				prob_left = 100 - base_prob
				new_prob = random.randint(0, prob_left)
				base_matrix.update({new_base:new_prob})
				base_prob = base_prob + new_prob
				bases.remove(new_base)
			new_base = random.choice(bases)
			prob_left = 100-base_prob
			base_matrix.update({new_base:prob_left})
			matrix.append(base_matrix)
		for base in mid_motif:
			matrix.append(RANDOM_SEQ)
		for base in end_motif:
			bases = ['A','T','C','G']
			base_matrix = {}
			base_prob = random.randint(base_conserv, 100)
			base_matrix.update({base:base_prob})
			bases.remove(base)
			for i in range(2):
				new_base = random.choice(bases)
				prob_left = 100 - base_prob
				new_prob = random.randint(0, prob_left)
				base_matrix.update({new_base:new_prob})
				base_prob = base_prob + new_prob
				bases.remove(new_base)
			new_base = random.choice(bases)
			prob_left = 100-base_prob
			base_matrix.update({new_base:prob_left})
			matrix.append(base_matrix)
	matrix = convertDistr(matrix)		
	return motif, matrix
	
# Create seeds for making footprints
def getSeeds(seed_num):
	seeds = []	
	motif_types = ['long', 'short', 'longcon', 'shortcon', 'split']	
	for i in range(seed_num):
		new_seed = []
		type = random.choice(motif_types)
		if type == 'long':	
			base_conserv, min_len, max_len, split = 55, 8, 16, 0
		elif type =='short':
			base_conserv, min_len, max_len, split = 55, 4, 8, 0
		elif type == 'longcon':
			base_conserv, min_len, max_len, split = 70, 8, 16, 0
		elif type == 'shortcon':
			base_conserv, min_len, max_len, split = 70, 4, 8, 0
		else:
			base_conserv, min_len, max_len, split = 70, 10, 16, 1
		motif, matrix = makeMotif(RANDOM_SEQ, base_conserv, min_len, max_len, split)
		new_seed.append(type)
		new_seed.append(motif)
		new_seed.append(matrix)
		seeds.append(new_seed)	
	return seeds	
	
# Create random nucleotides from a distribution	
def sampleDistr(distr):
	sample = random.randint(1,100)
	if sample <= distr.get('A'):
		return 'A'
	elif sample <= distr.get('T'):
		return 'T'
	elif sample <= distr.get('C'):
		return 'C'
	else:
		return 'G'

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
	
	
	

# main. Get info, generate seeds, make footprints, write footprints to output.
seed_count = int(sys.argv[1])
f_count = int(sys.argv[2])
foot_range = [int(sys.argv[3]),int(sys.argv[4])]
seeds = getSeeds(seed_count)
outfile = open(sys.argv[5], 'w')
for i in range(f_count):
	use_seed = random.choice(seeds)
	foot_length = random.randint(foot_range[0],foot_range[1])
	fprintWriter(use_seed[2], foot_length, use_seed[0], outfile, i, use_seed)
outfile.close()
	
	
	
	
