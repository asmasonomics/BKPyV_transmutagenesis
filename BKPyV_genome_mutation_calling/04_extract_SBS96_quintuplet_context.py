## python extract_SBS96_quintuplet_context.py file.tsv ref.fa

import sys
import pysam

tsv_in = open(sys.argv[1]).read().rstrip("\n").split("\n")
genome = pysam.FastaFile(sys.argv[2])

# take input TSV file, reverse string, split on '/' to get file name, reverse back again, split on first '.'
sampname = ((str(sys.argv[1])[::-1]).split("/")[0][::-1]).replace(".tsv", "")
outfile = open((sampname + "_SBS96_quintuplet_context.tsv"), "w")

# create list of chrom and position to check for neighboring sites
records = []
for i in tsv_in:
	records += [[i.split("\t")[0],i.split("\t")[1],i.split("\t")[2],i.split("\t")[3]]]

# loop through all sites
for record in records:
	# extract flanking bases to get quintuplet context
	flank_pre = genome.fetch(record[0], int(record[1])-3, int(record[1])-1)
	flank_post = genome.fetch(record[0], int(record[1]), int(record[1])+2)
	quin = (flank_pre + record[2] + flank_post).upper()

	# if mutation in poorly sequenced region, skip
	if "N" not in quin:
		# define alternate allele
		alt_allele = record[3]

		# COSMIC mut sigs all C/T, so reverse complement if G/A
		if record[2] in ["G", "A"]:
			quin = quin[::-1].lower().replace("a", "T").replace("c", "G").replace("t","A").replace("g","C")
			alt_allele = alt_allele.lower().replace("a", "T").replace("c", "G").replace("t","A").replace("g","C")

		# write original record and quin
		outfile.write(record[0] + "\t" + record[1] + "\t" + record[2] + "\t" + record[3] + "\t" + quin + "\t" + quin[1:4] + "\t" + quin[1] + "[" + quin[2] + ">" + alt_allele + "]" + quin[3] + "\n")

outfile.close()

