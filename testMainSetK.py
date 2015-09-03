# For the purpose of writing test scripts
# USE: python testing.py <input file location> <output file location>
# OUTPUT: Footprint sequence, sequence position(number), motif start position, Seed info (list)
import sys
import out
import main
import load

print ['test_start']
input_loc = str(sys.argv[1])
output_loc = str(sys.argv[2])
num_k = int(sys.argv[3])
out.writeOutput(output_loc, main.mainSetK(load.process_Generated(input_loc),num_k))
print ['test_end']
