import argparse
import ast
import math
import numpy as np
from collections import Counter
import json

#Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input file, need absolute path")
parser.add_argument("db_path", help="database directory")
#parser.add_argument("output_name", help="output file")
args = parser.parse_args()
input_file = args.input_file
db_path = args.db_path
#output = args.output_name

all_seq = []
fInput = open("%s"%input_file,'r')
for line in fInput:
	all_seq = all_seq + ast.literal_eval(line)
fInput.close()
all_seq = list(set(all_seq))

def mapping(db_file, mir_seqs):
	f = open("%s"%db_file,'r')
	db = []
	first_line = True
	for line in f:
		if line.startswith('>'):
			if first_line: first_line = False
			else: db.append(seq)
			seq = ""
		else:
			seq = seq.rstrip() + line.upper().replace("U","T").rstrip()
	db.append(seq)
	f.close()
	db = list(set(db))

	seq_match = []
	seq_unmatch = []
	for seq in mir_seqs:
		if seq in db:
			seq_match.append(seq)
		else:
			seq_unmatch.append(seq)

	return [seq_match, seq_unmatch]


if db_path.endswith('/'):
	db_file = db_path[:-1] + "/12fly-mir-flyb.fa"
else:
	db_file = db_path + "/12fly-mir-flyb.fa"
[db1_match, db1_unmatch] = mapping(db_file, all_seq)


if db_path.endswith('/'):
	db_file = db_path[:-1] + "/mature.fa"
else:
	db_file = db_path + "/mature.fa"
[db2_match, db2_unmatch] = mapping(db_file, db1_unmatch)

db_match = db1_match + db2_match

print db_match
