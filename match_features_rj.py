
import sys
import os
import re

WINDOW = 500

def open_and_parse_pairs(fp):
	pairs = []
	with open(fp) as f:
		for l in f: 
			# last two columns contains names of the peaks
			last_cols = l.strip().split('\t')[-2:]
			pairs.append(tuple(tuple(last_cols[0].split(', ')), tuple(last_cols[1].split(', '))))
	# [ ( ('peak_a1(1,2)', 'peak_a2(2,3)'), ('peak_b1(1,2)', 'peak_b2(2,3)') ), ( ... ), ... ]
	return pairs

# input: ('peak_a1(1,2)', 'peak_a2(2,3)'), ('peak_b1(1,2)', 'peak_b2(2,3)') )
# from these peaks, find the leftmost and rightmost boundaries
def get_min_max_window(peak_pairs):
	rx = r'[\(,](\d+)'
	left = min(re.findall(rx, peak_pairs[0][0])[0], re.findall(rx, peak_pairs[1][0])[0])
	right = max(re.findall(rx, peak_pairs[0][-1])[1], re.findall(rx, peak_pairs[1][-1])[1])
	return (left-WINDOW, right+WINDOW)

if __name__ == "__main__":
	del sys.argv[0]
	if not len(sys.argv) == 2:
		sys.stderr.write("Usage: [paired peaks tsv file] [features gff file]\n")
		sys.exit(1)
	else:
		paired_peaks_file = sys.argv[0]
		features_gff_file = sys.argv[1]
		if not (os.path.isfile(paired_peaks_file) and os.path.isfile(features_gff_file)):
			sys.stderr.write("Could not find file(s)")
			sys.exit(1)
		else:
			pairs = open_and_parse_pairs(paired_peaks_file)
			for pair in pairs:
				print(get_min_max_window(pair))