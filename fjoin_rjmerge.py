# fjoin_rjmerge.py
# ================
# takes fjoin output and chains together features that overlap; outputs the total feature 
# length added together along with the names of the features that are chained together
# 
# 2015.02.20 River Jiang
# 
# NOTE:
# make sure the COL_* constants reflect which columns the chr, start, stop, and names for each 
# paired feature can be found in the input data
# 
# fjoin outputs data sorted only by feature_start, and not by the chromosome
# make sure to sort data first by chromosome, and then by feature_start before feeding to this 
# script because this makes the assumption that feature chains are adjacent in the input data
# 	e.g.
# 	cat fjoin_condition.xls | sort -k2,2 -k3,3n > fjoin_condition_sorted.xls

import sys
import os

COL_COUNT = 9
COL_CHR = 0
COL_START = 1
COL_END = 2
COL_NAME = 8

def open_and_parse(fp):
	pairs = []
	with open(fp) as f:
		for l in f:
			ls = l.strip().split("\t")
			pairs.append([ ( ls[1+COL_CHR], ls[1+COL_START], ls[1+COL_END], ls[1+COL_NAME] ), \
				( ls[1+COL_CHR+COL_COUNT], ls[1+COL_START+COL_COUNT], ls[1+COL_END+COL_COUNT], ls[1+COL_NAME+COL_COUNT] ) ])
	return pairs

def get_index(s):
	return s[3].split("_")[-1]

def get_linked_pairs(pairs):
	linked_pairs = []
	curr_links = [ pairs[0] ]

	for curr_pair in range(len(pairs) - 1):
		search_pair = curr_pair+1
		if (get_index(pairs[search_pair][0]) > get_index(pairs[curr_pair][0])) and \
		   (get_index(pairs[search_pair][1]) > get_index(pairs[curr_pair][1])):
			linked_pairs.append(curr_links)
			curr_links = [ pairs[curr_pair + 1] ]
		else:
			curr_links.append(pairs[search_pair])

	linked_pairs.append(curr_links)
	return linked_pairs

def get_unique_sets(linked_pairs):
	unique_sets = []
	for chain in linked_pairs:
		xset = set()
		yset = set()
		for pair in chain:
			xset.add(pair[0])
			yset.add(pair[1])
		unique_sets.append([ list(xset), list(yset) ])
	return unique_sets

def make_summary(unique_sets):
	summary = ""
	for chain in unique_sets:
		sorted_chain = [ sorted(chain[0], key=lambda x: x[1]), sorted(chain[1], key=lambda x: x[1]) ]

		xsize = 0
		xnames = []
		for x in sorted_chain[0]:
			xsize += (int(x[2]) - int(x[1]))
			xnames.append(x[3] + "(" + x[1] + "," + x[2] + ")")

		ysize = 0
		ynames = []
		for y in sorted_chain[1]:
			ysize += (int(y[2]) - int(y[1]))
			ynames.append(y[3] + "(" + y[1] + "," + y[2] + ")")
		
		summary += chain[0][0][0] + "\t" + str(xsize) + "\t" + str(ysize) + "\t" + ", ".join(xnames) + "\t" + ", ".join(ynames) + "\n"
	return summary

if __name__ == "__main__":
	del sys.argv[0]
	if len(sys.argv) != 1:
		sys.stderr.write("no file given\n")
		sys.exit(1)
	else:
		fp = sys.argv[0]
		if not os.path.isfile(fp):
			sys.stderr.write("could not find file\n")
			sys.exit(1)
		else:
			pairs = open_and_parse(fp)
			linked_pairs = get_linked_pairs(pairs)
			unique_sets = get_unique_sets(linked_pairs)
			print(make_summary(unique_sets))