import sys
import os
import re

fastq_regex = r'Meneghini_[0-9+]_(.+)_([ATGC]{6})_L001_R1_001.fastq.gz'

RGSM = "KY_ChIP"
RGCN = "CCBR"

bowtie2indexpath = "../Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/Bowtie2Index/genome"
picardtoolspath = "../picard-tools-1.119/"
gatkpath = "../GenomeAnalysisTK.jar"
gatkgenome = "../sacCer3/genome.fa"

# check needed files are present
for f in [os.path.dirname(bowtie2indexpath), picardtoolspath, gatkpath, gatkgenome]:
	if not os.path.exists(f):
		print("Check config: " + f + " does not exist")
		exit(1)

if len(sys.argv) < 2:
	print("""Usage: sam_to_bam.py [folder] -p [prefix]

Runs MX sam to bam pipeline. Feed me path to folder containing fastq.gz files, and
I will output bam files into the current directory.""")
	exit(1)

del(sys.argv[0])
path = sys.argv[0]
del(sys.argv[0])
if not os.path.exists(path):
	print("Folder does not exist.")
	exit(1)

prefix = ""
if len(sys.argv) == 2:
	if sys.argv[0] == "-p":
		del sys.argv[0]
		prefix = sys.argv[0] + "_"

fastq_files = []
for file in os.listdir(path):
	if file.endswith(".fastq.gz"):
		fastq_files.append(file)

if len(fastq_files) == 0:
	print("No files found. Check files exist, and check fastq_regex in config.")
	exit(1)

prettynames = []
rgpus = []
for f in fastq_files:
	m = re.match(fastq_regex, f)
	prettynames.append(m.group(1))
	rgpus.append(m.group(2))

for i in range(len(fastq_files)):
	print(path + "/" + fastq_files[i] + ": " + prettynames[i] + " -> " + rgpus[i])

raw_input("\nContinue...?")

commands = ""

# bowtie2 map to genome
for i in range(len(fastq_files)):
	commands += "bowtie2 -x " + bowtie2indexpath + " -U " + path + "/" + fastq_files[i] + " -S " + prefix + prettynames[i] + ".sam -p 4\n"

# convert sam to bam
for i in range(len(fastq_files)):
	commands += "java -Xmx2g -jar " + picardtoolspath + "SamFormatConverter.jar INPUT=" + prefix + prettynames[i] + ".sam OUTPUT=" + prefix + prettynames[i] + ".bam\n"

# add or replace read groups
for i in range(len(fastq_files)):
	commands += "java -Xmx2g -jar " + picardtoolspath + "AddOrReplaceReadGroups.jar INPUT=" + prefix + prettynames[i] + ".bam OUTPUT=" + prefix + prettynames[i] + "_ordered.bam SORT_ORDER=coordinate RGPL=illumina RGPU=" + rgpus[i] + " RGSM=" + RGSM + " RGCN=" + RGCN + " RGLB=NA\n"

# index bam files
for i in range(len(fastq_files)):
	commands += "java -Xmx2g -jar " + picardtoolspath + "BuildBamIndex.jar INPUT=" + prefix + prettynames[i] + "_ordered.bam\n"

# depth of coverage
for i in range(len(fastq_files)):
	commands += "java -Xmx2g -jar " + gatkpath + " -R "+ gatkgenome + " -T DepthOfCoverage -o " + prefix + prettynames[i] + "_DoC.sgr -I " + prefix + prettynames[i] + "_ordered.bam --omitIntervalStatistics --omitLocusTable\n"

print commands

os.system(commands)