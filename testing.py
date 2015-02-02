#For the purpose of writing test scripts

outputLoc = 'testData/generated1'

def generate_standard():
	os.system(("python motifMaker.py 100 1000 5 14 " + outputLoc))
	return

def test_generated():
	data = Load.process_Generated()
	peaks = get_seqs (data,chrom)
	for peak in peaks:
		print peak
	return

print ['test_start']

print ['test_end']
