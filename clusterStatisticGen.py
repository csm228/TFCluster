
# A script for generating statistics on clustered generated data
import sys
import math
import random

numRandReplicates = 1000

#loads the cluster data from a single cluster file
def loadClusterData (clusterFile):
	source = open(clusterFile)
	clusters = []
	totalDict = {}
	sourceStrings = source.read().split('Cluster')
	#need to use sourceStrings[1:] b/c split adds an empty string before the first entry
	for potentialCluster in sourceStrings[1:]:
		clusterLines = potentialCluster.split('\n')
		IDname = clusterLines [0]
		description = clusterLines[1]
		description = clusterLines[2]
		memberMotifs = []
		numMembers = 0
		for memberSeq in clusterLines[3:]:
			elementData = memberSeq.split('\t')
			if len(elementData) > 1:
				motif = elementData[1]
				# if motif == '':
				# 	print 'empty motif'
				# print motif
				memberMotifs += [motif]
				numMembers += 1
				if motif in totalDict:
					totalDict[motif] = totalDict[motif] + 1
				else:
					totalDict[motif] = 1
		if numMembers > 0:
			clusters += [(numMembers,memberMotifs)]
	return (clusters, totalDict)

def analyze (clusterData):
	(clusters, totalDict) = clusterData
	homogeneities = []
	for (numMembers,memberMotifs) in clusters:
		clusterMotifDict = {}
		for motif in memberMotifs:
			if motif in clusterMotifDict:
				clusterMotifDict[motif] = clusterMotifDict[motif] + 1
			else:
				clusterMotifDict[motif] = 1
		motifFreq = clusterMotifDict.items()
		motifFreq.sort(key = lambda seg: seg[1], reverse=True)
		# print motifFreq[0]
		clusterHomogeneity = float(motifFreq[0][1])/float(numMembers)
		homogeneities += [clusterHomogeneity]
	homogeneities.sort(reverse=True)
	return homogeneities

# def predictNull ():
# 	return

#excel format? rows = individual clusterings
def writeHomogeneitiesListToTable(homogeneitiesList, outfileName):
	outfile = open(outfileName, 'w')
	clusterNum = 1
	for homogeneities in homogeneitiesList:
		line = "Dataset" + str(clusterNum)
		for percent in homogeneities:
			line += "\t" + str(percent)
		line += "\n"
		outfile.write(line)
		clusterNum += 1
	outfile.close()
	return

def categorizedHomogeneities(homogeneitiesList):
	numG75 = 0
	numG50 = 0
	numL50 = 0
	total = 0
	for homogeneities in homogeneitiesList:
		for percent in homogeneities:
			if percent > 0.75:
				print 'G75'
				numG75 += 1
			elif percent > 0.50:
				print 'G50'
				numG50 += 1
			else:
				print 'L50'
				numL50 += 1
			total += 1
	portionG75 = float(numG75)/float(total)
	portionG50 = float(numG50)/float(total)
	portionL50 = float(numL50)/float(total)
	return (portionG75,portionG50,portionL50)

#test null by resampling clusters of the given size

def generateRandomHomogeneities(clusterData):
	(clusters, totalDict) = clusterData
	homogeneities = []
	allMotifs = []
	for i in range(numRandReplicates):
		for (numMembers,memberMotifs) in clusters:
			allMotifs += memberMotifs
		random.shuffle(allMotifs)
		position = 0
		for (numMembers,memberMotifs) in clusters:
			# print numMembers
			# print position
			clusterMotifDict = {}
			for motif in allMotifs[position:(position+numMembers)]:
				if motif in clusterMotifDict:
					clusterMotifDict[motif] = clusterMotifDict[motif] + 1
				else:
					clusterMotifDict[motif] = 1
			# print clusterMotifDict
			motifFreq = clusterMotifDict.items()
			motifFreq.sort(key = lambda seg: seg[1], reverse=True)
			# print motifFreq
			clusterHomogeneity = float(motifFreq[0][1])/float(numMembers)
			homogeneities += [clusterHomogeneity]
			position = position + numMembers
	homogeneities.sort(reverse=True)
	return homogeneities

#writes an R table of data and randomized homogeneities
def writeHomogDistribTable(homogeneities,randHomogeneities,outfileName):
	outfile = open(outfileName, 'w')
	outfile.write("Clustered,Randomized\n")
	numHom = len(homogeneities)
	for i in range(len(randHomogeneities)):
		# % numRandReplicates to account for the replicate sampling in the randomized data
		outfile.write(str(homogeneities[int(i % numHom)])+','+str(randHomogeneities[i])+'\n')
	outfile.close()
	return

#Parameters: file location/name, number of datasets, number of runs
dataName = str(sys.argv[1])
numDataSets = int(sys.argv[2])
numRuns = int(sys.argv[3])
outfileName = str(sys.argv[4])




homogeneitiesList = []
randomHomogeneitiesList = []
homogeneities = []
randHomogeneities = []
for i in range(numDataSets):
	for j in range(numRuns):
		clusterData = loadClusterData(dataName + str(i+1) + 'v' + str(j+1))
		# homogeneitiesList += [analyze(clusterData)]
		homogeneities += analyze(clusterData)
		# randomHomogeneitiesList += [generateRandomHomogeneities(clusterData)]
		randHomogeneities += generateRandomHomogeneities(clusterData)
# writeHomogeneitiesListToTable(homogeneitiesList,outfileName)
# print categorizedHomogeneities(homogeneitiesList)
writeHomogDistribTable(homogeneities,randHomogeneities,outfileName)
		

# print(analyze(loadClusterData(dataName + '1v1')))








