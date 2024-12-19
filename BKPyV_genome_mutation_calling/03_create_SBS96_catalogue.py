## python create_SBS96_catalogue.py file.tsv ref.fa

import sys
import pysam

vcf_in = open(sys.argv[1]).read().rstrip("\n").split("\n")
genome = pysam.FastaFile(sys.argv[2])

# take input tsv, reverse string, split on '/' to get file name, reverse back again, split on first '.'
sampname = ((str(sys.argv[1])[::-1]).split("/")[0][::-1]).replace(".tsv", "")
outfile = open((sampname + "_SBS96_catalogue.tsv"), "w")

# loop through records to process SNPs
poss_muts = [("CA","ACA"), ("CA","ACC"), ("CA","ACG"), ("CA","ACT"), ("CA","CCA"), ("CA","CCC"), ("CA","CCG"), ("CA","CCT"), ("CA","GCA"), ("CA","GCC"), ("CA","GCG"), ("CA","GCT"), ("CA","TCA"), ("CA","TCC"), ("CA","TCG"), ("CA","TCT"), ("CG","ACA"), ("CG","ACC"), ("CG","ACG"), ("CG","ACT"), ("CG","CCA"), ("CG","CCC"), ("CG","CCG"), ("CG","CCT"), ("CG","GCA"), ("CG","GCC"), ("CG","GCG"), ("CG","GCT"), ("CG","TCA"), ("CG","TCC"), ("CG","TCG"), ("CG","TCT"), ("CT","ACA"), ("CT","ACC"), ("CT","ACG"), ("CT","ACT"), ("CT","CCA"), ("CT","CCC"), ("CT","CCG"), ("CT","CCT"), ("CT","GCA"), ("CT","GCC"), ("CT","GCG"), ("CT","GCT"), ("CT","TCA"), ("CT","TCC"), ("CT","TCG"), ("CT","TCT"), ("TA","ATA"), ("TA","ATC"), ("TA","ATG"), ("TA","ATT"), ("TA","CTA"), ("TA","CTC"), ("TA","CTG"), ("TA","CTT"), ("TA","GTA"), ("TA","GTC"), ("TA","GTG"), ("TA","GTT"), ("TA","TTA"), ("TA","TTC"), ("TA","TTG"), ("TA","TTT"), ("TC","ATA"), ("TC","ATC"), ("TC","ATG"), ("TC","ATT"), ("TC","CTA"), ("TC","CTC"), ("TC","CTG"), ("TC","CTT"), ("TC","GTA"), ("TC","GTC"), ("TC","GTG"), ("TC","GTT"), ("TC","TTA"), ("TC","TTC"), ("TC","TTG"), ("TC","TTT"), ("TG","ATA"), ("TG","ATC"), ("TG","ATG"), ("TG","ATT"), ("TG","CTA"), ("TG","CTC"), ("TG","CTG"), ("TG","CTT"), ("TG","GTA"), ("TG","GTC"), ("TG","GTG"), ("TG","GTT"), ("TG","TTA"), ("TG","TTC"), ("TG","TTG"), ("TG","TTT")]
CA_trip = []
CG_trip = []
CT_trip = []
TA_trip = []
TC_trip = []
TG_trip = []

# create list of chrom and position to check for neighboring sites
vcf_records = []
for i in vcf_in:
	vcf_records += [[i.split("\t")[0],i.split("\t")[1],i.split("\t")[2],i.split("\t")[3]]]
vcf_records_len = len(vcf_records)-1

numrecords = 0
for record in vcf_records:
	
	# check if SNVs are neighbouring. If so, don't include in analysis.	
	include = True	
	if numrecords == 0:
		if ((vcf_records[0][0] == vcf_records[1][0]) and (int(vcf_records[0][1])+1 == int(vcf_records[1][1]))):
			include = False

	elif numrecords == vcf_records_len:
		if ((vcf_records[numrecords-1][0] == vcf_records[numrecords][0]) and (int(vcf_records[numrecords-1][1])+1 == int(vcf_records[numrecords][1]))):
			include = False

	else:
		if (((vcf_records[numrecords-1][0] == vcf_records[numrecords][0]) and (int(vcf_records[numrecords-1][1])+1 == int(vcf_records[numrecords][1]))) or ((vcf_records[numrecords][0] == vcf_records[numrecords+1][0]) and (int(vcf_records[numrecords][1])+1 == int(vcf_records[numrecords+1][1])))):
			include = False


	numrecords += 1

	# uncomment the line below to include all sites, regardless of whether SNVs are next to each other
	#include = True

	# neighbouring sites will change include to False and won't be included from here
	if include:
		flank_pre = genome.fetch(record[0], int(record[1])-2, int(record[1])-1)
		flank_post = genome.fetch(record[0], int(record[1]), int(record[1])+1)

		triplet = (flank_pre + record[2] + flank_post).upper()
		if "N" not in triplet:
			alt_allele = record[3]
			if record[2] in ["G", "A"]:
				triplet = triplet[::-1].lower().replace("a", "T").replace("c", "G").replace("t","A").replace("g","C")
				alt_allele = alt_allele.lower().replace("a", "T").replace("c", "G").replace("t","A").replace("g","C")
			if triplet[1] == "C":
				if alt_allele == "A":
					CA_trip += [triplet]
				elif alt_allele == "G":
					CG_trip += [triplet]
				elif alt_allele == "T":
					CT_trip += [triplet]
			elif triplet[1] == "T":
				if alt_allele == "A":
					TA_trip += [triplet]
				elif alt_allele == "C":
					TC_trip += [triplet]
				elif alt_allele == "G":
					TG_trip += [triplet]

	

from itertools import groupby
def list_sorter_counter(my_list):
	my_list.sort()
	count = [(key, len(list(group))) for key, group in groupby(my_list)]
	return count

counts = {}
for i,k in [(CA_trip,"CA"), (CG_trip,"CG"), (CT_trip,"CT"), (TA_trip,"TA"), (TC_trip,"TC"), (TG_trip,"TG")]:
	for j in list_sorter_counter(i):
		#print(j[0],round(((int(j[1])*100)/numrecords),2))
		#print(str(j[0]) + "\t" + str(j[1]))
		counts[(k, str(j[0]))] = str(j[1])

outfile.write("Substitution\t" + sampname + "\n")
for mut in poss_muts:
	if mut in counts:
		#print(mut, counts[mut])
		outfile.write(mut[1][0] + "[" + mut[0][0]+ ">" + mut[0][1] + "]" + mut[1][-1] + "\t" + str(counts[mut]) + "\n")
	else:
		#print(mut, 0)
		outfile.write(mut[1][0] + "[" + mut[0][0]+ ">" + mut[0][1] + "]" + mut[1][-1] + "\t0\n")
outfile.close()


