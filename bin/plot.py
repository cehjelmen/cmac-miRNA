import argparse
import ast
import math
import matplotlib.pylab as plt
import numpy as np
from collections import Counter
import json

#Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("input_file", help="input file, need absolute path")
parser.add_argument("output_name", help="output file")
args = parser.parse_args()
input_file = args.input_file
output_name = args.output_name

fInput = open("%s"%input_file,'r')

#fInput = open("sample_trim.json",'r')
line = fInput.readline()
seq_dict = ast.literal_eval(line)
fInput.close()

unique_len_freq = dict()
for leng in seq_dict:
	unique_len_freq[int(leng)] = len(seq_dict[leng])

lists = sorted(unique_len_freq.items()) 
x, y = zip(*lists) 
plt.bar(x, y)
plt.xlabel("sequence length")
plt.ylabel("number of unique sequence")
plt.title("length distribution of unique sequence")
plt.savefig(output_name + "_unique_len_freq")
plt.clf()

all_len_freq = dict()
for leng in seq_dict:
	all_len_freq[int(leng)] = 0
	for seq in seq_dict[leng]:
		all_len_freq[int(leng)] += seq_dict[leng][seq]

lists = sorted(all_len_freq.items()) 
x, y = zip(*lists) 
plt.bar(x, y)
plt.xlabel("sequence length")
plt.ylabel("number of all sequence")
plt.title("length distribution of all sequence")
plt.savefig(output_name + "_all_len_freq")
plt.clf()