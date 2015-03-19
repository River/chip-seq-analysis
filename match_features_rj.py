# match_features_rj.py
# ====================
# takes peak pairs output from fjoin_rjmerge.py as well as annotated features
# gff file, and outputs all features that are within +/- WINDOW size of the
# leftmost and rightmost edges of the broadest peak of each pair
# 
# can manually combine output with the paired peaks in a big excel spreadsheet
# for GO enrichment analysis, etc.

import sys
import os
import re

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
			# ignore empty lines
			if len(l) > 1:
				# last two columns contains names of the peaks
				cols = l.strip().split('\t')
				pairs.append(
					tuple(
						[
							cols[0], 
							tuple(cols[3].split(', ')), 
							tuple(cols[4].split(', ')), 
							cols[1], 
							cols[2]
						]
					)
				)
	# [ ( 'chr1', ('peak_a1(1,2)', 'peak_a2(2,3)'), ('peak_b1(1,2)', 'peak_b2(2,3)'), widthA, widthB ), ( ... ), ... ]
	return pairs

# input: ( 'chr1', ('peak_a1(1,2)', 'peak_a2(2,3)'), ('peak_b1(1,2)', 'peak_b2(2,3)') )
# from these peaks, find the leftmost and rightmost boundaries
# assume that peaks are sorted!! will only look at leftmost and rightmost elems in tuples
def get_min_max_window(peak_pairs):
	# match numbers that come after a ( or a ,
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
def get_features_in_window(chrom, window, overlap_mode):
	# overlap_mode:
	#   0 start (TSS)
	#	1 stop (TTS)
	#	2 both start and stop
	features = []
	for f in FEATURES:
		# check if chromosome is correct
		if f[0] == chrom:
			# window overlaps TSS
			if overlap_mode == 0:
				if (f[1] >= window[0]) and (f[1] <= window[1]):
					features.append(f[3])

			# window overlaps TTS
			elif overlap_mode == 1:
				if (f[2] >= window[0]) and (f[2] <= window[1]):
					features.append(f[3])
			
			# window overlaps both TSS and TTS
			elif overlap_mode == 2:
				# overlap: min(endA, endB) - max(startA, startB) > 0
				if (min(window[1], f[2]) - max(window[0], f[1])) > 0:
					features.append(f[3])

			# if we have gone past the window, stop looking
			elif f[1] >= window[1]:
				break
	return features

if __name__ == "__main__":
	del sys.argv[0]
	if not len(sys.argv) == 4:
		sys.stderr.write("Usage: [overlap_window] [overlap_mode] [paired peaks tsv file from fjoin_rjmerge.py] [features gff file]\n")
		sys.exit(1)
	else:
		WINDOW = int(sys.argv[0])
		del sys.argv[0]
		sys.stderr.write("Looking for features +/- " + str(WINDOW) + " of peaks\n")

		overlap_mode = int(sys.argv[0])
		del sys.argv[0]
		foo = ""
		if overlap_mode == 0:
			foo = "TSS"
		elif overlap_mode == 1:
			foo = "TTS"
		elif overlap_mode == 2:
			foo = "TSS or TTS"
		sys.stderr.write("Looking for overlap of the " + foo + "\n")

		paired_peaks_file = sys.argv[0]
		features_gff_file = sys.argv[1]
		if not (os.path.isfile(paired_peaks_file) and os.path.isfile(features_gff_file)):
			sys.stderr.write("Could not find file(s)\n")
			sys.exit(1)
		else:
			FEATURES = open_and_parse_gff(features_gff_file)

			# OUTPUT
			# chr,totalwidthA,totalwidthB,peaksA,peaksB,window_size,features
			summary = ""
			pairs = open_and_parse_pairs(paired_peaks_file)
			for pair in pairs:
				chrom = pair[0]
				window = get_min_max_window(pair)
				features = get_features_in_window(chrom, window, overlap_mode)
				for f in features:
					summary += chrom + "\t" + str(WINDOW) + "\t" + pair[3] + "\t" + pair[4] + "\t" + ", ".join(pair[1]) + "\t" + ", ".join(pair[2]) + "\t" + f + "\n"

			print(summary)
