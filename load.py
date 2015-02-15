import ast

#For the comparison during alignment, ensure uppercase in peak strings

#takes a generated test data file and outputs a list of peaks (lists)
def process_Generated (testData):
	source = open(testData)
	data = []
	for line in source:
		datum = (line.strip()).split('\t')
		datum[0] = (str(datum[0])).upper()
		datum[1] = int(datum[1])
		datum[2] = int(datum[2])
		datum[3] = ast.literal_eval(datum[3])
		# datum[3] = [e.strip() for e in datum[3]]
		#data.append(datum)
		data += [datum]
	#print data[0]
	#print data[4]
	return data

#takes an ENCODE narrowPeak file and outputs a list of lists (peaks)
def process_DNase (narrowPeak):
	source = open(narrowPeak)
	data = []
	for line in source:
		datum = line.strip().split('\t')
		datum[1] = int(datum[1])
		datum[2] = int(datum[2])
		datum[4] = int(datum[4])
		datum[6] = int(datum[6])
		datum[7] = float(datum[7])
		datum[8] = float(datum[8])
		datum[9] = float(datum[9])
		data.append(datum)
	return data

#takes an ENCODE chromosome file and outputs a string of its sequence
def process_chrom (chromosome):
	raw = open(chromosome)
	raw.readline()
	seq = ''
	for line in raw:
		seq = seq + line.strip()
	return seq


def get_datum(data,i):
	return data[i]

def get_chrom(datum):
	return datum[1]

def get_sequence(datum, seq):
	return seq[datum[1]:datum[2]]

def get_seqs(data,chromosome):
	peaks = []
	for datum in data:
		peaks.append(get_sequence(datum,chromosome))
	return peaks


def test():
	chrom = process_chrom('data/chr21.fa')
	data = process_DNase('data/ENCFF001WIR.narrowPeak_chr21.np')
	peaks = get_seqs (data,chrom)
	for peak in peaks:
		print peak
	return


print ['load_test - odd python attributes']
#test()
print ['load_test2']