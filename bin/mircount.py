import argparse
import ast
import math
import numpy as np
from collections import Counter
import json

#Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input file, need absolute path")
parser.add_argument("mir_path", help="database directory")
#parser.add_argument("output_name", help="output file")
args = parser.parse_args()
input_file = args.input_file
mir_path = args.mir_path
#output = args.output_name

fInput = open("%s"%input_file,'r')
line = fInput.readline()
seq_dict = ast.literal_eval(line)
fInput.close()

mir_seq = []
fInput = open("%s"%mir_path,'r')
line = fInput.readline()
mir_seq = ast.literal_eval(line)
fInput.close()

for seq in mir_seq:
	if seq in seq_dict[str(len(seq))]:
		print seq, ",", seq_dict[str(len(seq))][seq]
	else:
		print seq, ",", 0