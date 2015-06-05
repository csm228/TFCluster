# For the purpose of writing test scripts
# USE: python testing.py <input file location> <output file location>
# OUTPUT: Footprint sequence, sequence position(number), motif start position, Seed info (list)
import sys
import load
import out
import main

import align
import seed
import paring
import cluster

import matrixSeed

peak1 = ['ATTGCTCTAGCTACGATCTATACTGTACACAAA',0]
peak2 = ['GCGCGTAATCATTACGTATTTCCGGCGATATCGGGGGG',1]
peak3 = ['GCGCGTAGGCGCTATTTCTG',2]
peak4 = ['TTCTTTTTTCGCTTATATCGGCTAT',3]
peak5 = ['TCCTCTACTG',4]
peak6 = ['CGTATTATCTTTCGTTTGTCTTATCGTTTCTTATTTCGCGCGATTCTAGCTTATGCTATTGCTA',5]

mean1 = []
meanWords1 = []

def test_seed_abstraction1():
	print 'seed: ' + str(peak1)
	mean1 = seed.abstract(peak1)
	print 'abstracted seed: ' + str(mean1)
	return mean1

def test_wordification1():
	print 'mean: ' + str(mean1)
	(mean,targetLength) = mean1
	meanWords1 = align.wordify(mean)
	print 'mean words: ' + str(meanWords1)
	return meanWords1

def test_align(peak1, peak2):
	(mean, targetLength) = seed.abstract(peak2)
	print align.align(peak1,align.wordify(mean),targetLength)

def test_align1():
	print 'peak: ' + str(peak1)
	print 'mean: ' + str(mean1)
	(i, score) = align.align(peak1,meanWords1,-1)
	print 'index is ' + str(i) + ', should be 0?'
	print 'score is ' + str(score)
	return

def test_pickInitMeans():
	peaks = [peak1,peak2,peak3,peak4,peak5,peak6]
	print seed.pickInitMeans(peaks,3)

def test_pickInitSeeds():
	peaks = [peak1,peak2,peak3,peak4,peak5,peak6]
	print matrixSeed.pickInitSeeds(peaks,3)


print 'test_start \n'
# mean1 = test_seed_abstraction1()
print '\n'
# meanWords1 = test_wordification1()
print '\n'
# test_align(peak1,peak1)
# test_align(peak2,peak2)
# test_align(peak3,peak3)
print '\n'
# test_pickInitMeans()
print '\n'
test_pickInitSeeds()
print '\n'
# input_loc = str(sys.argv[1])
# output_loc = str(sys.argv[2])
# out.writeOutput(output_loc, main.main(load.process_Generated(input_loc)))
print 'test_end'
