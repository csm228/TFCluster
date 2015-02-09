# For the purpose of writing test scripts
# USE: python testing.py <input file location> <output file location>
# OUTPUT: Footprint sequence, sequence position(number), motif start position, Seed info (list)
import sys
import out
import main
import load

def test_generated(filename):
	data = Load.process_Generated()
	peaks = get_seqs (data,chrom)
	for peak in peaks:
		print peak
	return

print ['test_start']
input_loc = str(sys.argv[1])
output_loc = str(sys.argv[2])
out.writeOutput(output_loc, main.main(load.process_Generated(input_loc)))
print ['test_end']
