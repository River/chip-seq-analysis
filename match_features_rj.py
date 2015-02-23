
import sys
import os
import re

# how many bp upstream and downstream to look for features
WINDOW = 500

# don't edit below!
chr_conversion = {
	'chr1': 	'chrI',
	'chr2': 	'chrII',
	'chr3': 	'chrIII',
	'chr4': 	'chrIV',
	'chr5': 	'chrV',
	'chr6': 	'chrVI',
	'chr7': 	'chrVII',
	'chr8': 	'chrVIII',
	'chr9': 	'chrIX',
	'chr10': 	'chrX',
	'chr11': 	'chrXI',
	'chr12': 	'chrXII',
	'chr13': 	'chrXIII',
	'chr14': 	'chrXIV',
	'chr15': 	'chrXV',
	'chr16': 	'chrXVI',
	'chrM': 	'chrM'
}
FEATURES = ()

def open_and_parse_pairs(fp):
	pairs = []
	with open(fp) as f:
		for l in f: 
			# last two columns contains names of the peaks
			chrom = l.split('\t')[0]
			last_cols = l.strip().split('\t')[-2:]
			pairs.append(tuple([chrom, tuple(last_cols[0].split(', ')), tuple(last_cols[1].split(', '))]))
	# [ ( 'chr1', ('peak_a1(1,2)', 'peak_a2(2,3)'), ('peak_b1(1,2)', 'peak_b2(2,3)') ), ( ... ), ... ]
	return pairs

# input: ( 'chr1', ('peak_a1(1,2)', 'peak_a2(2,3)'), ('peak_b1(1,2)', 'peak_b2(2,3)') )
# from these peaks, find the leftmost and rightmost boundaries
# assume that peaks are sorted!! will only look at leftmost and rightmost elems in tuples
def get_min_max_window(peak_pairs):
	rx = r'[\(,](\d+)'
	left = min(int(re.findall(rx, peak_pairs[1][0])[0]), int(re.findall(rx, peak_pairs[2][0])[0]))
	right = max(int(re.findall(rx, peak_pairs[1][-1])[1]), int(re.findall(rx, peak_pairs[2][-1])[1]))
	return (left-WINDOW, right+WINDOW)

# FEATURES tuple:
# ( ( chr, start, end, feature_name ), ( ... ) , ... )
def open_and_parse_gff(fp):
	features = []
	with open(fp) as f:
		for l in f:
			if not l[0] == '#':
				cols = l.strip().split('\t')
				features.append(tuple([chr_conversion[cols[0]], int(cols[3]), int(cols[4]), cols[8]]))
	return tuple(features)

# return all features for which there is any overlap with the window
# assume that FEATURES contains all features, and that these are SORTED!
# will stop looking once it reaches a point where a feature is completely to the right of the window
def get_features_in_window(chrom, window):
	features = []
	for f in FEATURES:
		if f[0] == chrom:
			# overlap: min(endA, endB) - max(startA, startB) > 0
			if (min(window[1], f[2]) - max(window[0], f[1])) > 0:
				features.append(f[3])
			# if we have gone past the window, stop looking
			if f[1] >= window[1]:
				break
	return features

if __name__ == "__main__":
	del sys.argv[0]
	if not len(sys.argv) == 2:
		sys.stderr.write("Usage: [paired peaks tsv file] [features gff file]\n")
		sys.exit(1)
	else:
		paired_peaks_file = sys.argv[0]
		features_gff_file = sys.argv[1]
		if not (os.path.isfile(paired_peaks_file) and os.path.isfile(features_gff_file)):
			sys.stderr.write("Could not find file(s)\n")
			sys.exit(1)
		else:
			FEATURES = open_and_parse_gff(features_gff_file)

			pairs = open_and_parse_pairs(paired_peaks_file)
			for pair in pairs:
				chrom = pair[0]
				window = get_min_max_window(pair)
				features = get_features_in_window(chrom, window)
				print(pair)
				print(window)
				print(features)