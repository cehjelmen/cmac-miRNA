import argparse
import ast
import math
import numpy as np
from collections import Counter
import json

#Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input file, need absolute path")
parser.add_argument("min_len", help="min mirna seq length")
parser.add_argument("max_len", help="max mirna seq length")
#parser.add_argument("output_name", help="output file")
args = parser.parse_args()
input_file = args.input_file
min_len = args.min_len
max_len = args.max_len
#output = args.output_name

fInput = open("%s"%input_file,'r')
#fInput = open("sample_trim.json",'r')
line = fInput.readline()
seq_dict = ast.literal_eval(line)
fInput.close()

mir_list = []
rna_list = []
for length in seq_dict:
	if int(length) < int(min_len) or int(length) > int(max_len):
		rna_list = rna_list + seq_dict[length].keys()
	else:
		mir_list = mir_list + seq_dict[length].keys()

print mir_list