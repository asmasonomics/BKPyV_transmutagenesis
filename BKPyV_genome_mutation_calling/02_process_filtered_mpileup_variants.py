# python process_filtered_mpileup_variants.py sample.MD.REsitefilt.noindels.accepted-variants.tsv

import sys

chr = ""
pos = ""
freq = 0
ref = ""
alt = ""

# create output file
out = open(((sys.argv[1]).replace("tsv", "chr-pos-ref-alt.tsv")), "w")

# read the filtered var files
vars = open(sys.argv[1]).read().rstrip("\n").split("\n")

for site in vars:
    ## the replace here is for the BKPyV Dunlop genome (the "." in ".1" is lost in the mpileup processing)
	chr = (site.split("\t")[0]).replace("NC_0015381", "NC_001538.1")
	pos = site.split("\t")[1]
	ref = site.split("\t")[2]
	alt = site.split("\t")[-1]
	freq = float(site.split("\t")[4])

    # frequency check is to avoid carry through of non-reference assembly variants in viral population
	if freq <= 0.05:
		if len(alt) > 1:
			for nonref in alt:
				out.write(chr + "\t" + pos + "\t" + ref + "\t" + nonref + "\n")

		else:
			out.write(chr + "\t" + pos + "\t" + ref + "\t" + alt + "\n")

out.close()
