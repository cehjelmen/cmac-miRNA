import argparse
import numpy as np
from collections import Counter
import json

#Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input file, need absolute path")
parser.add_argument("output_file", help="output file")
args = parser.parse_args()
input_file = args.input_file
output_file = args.output_file

seq_list = list()
seq_dict = dict()
seq_uniq = dict()

fInput = open("%s"%input_file,'r')
flag = 0;
for line in fInput:
	flag += 1
	if (flag % 4 == 2):
		if line.endswith('\n'): line = line[:-1]
		line = line.replace("N", "")
		seq_list.append(line)
fInput.close()

for seq in seq_list:
	if len(seq) not in seq_dict:
		seq_dict[len(seq)] = list()
	seq_dict[len(seq)].append(seq)

for length in seq_dict:
	seq_uniq[length] = (dict(Counter(seq_dict[length])))

output_name = output_file + ".json"
with open("%s"%output_name,'w') as f:
	json.dump(seq_uniq, f)
