# python 05_APOBEC3A-APOBEC3B_quintuplet_summary.py /dirpath/*_SBS96_quintuplet_context.tsv

# import libraries
import sys
import os
import pandas as pd
from collections import Counter

# generate list of files
quin_files = [x for x in os.listdir(sys.argv[1]) if x.endswith("_SBS96_quintuplet_context.tsv")]
quin_files.sort()

# loop through datafiles, storing row outputs to make dataframe
qf_out = []
apobec_pent_types = ["RT[C>A]A", "RT[C>A]T", "RT[C>G]A", "RT[C>G]T", "RT[C>T]A", "RT[C>T]T", "YT[C>A]A", "YT[C>A]T", "YT[C>G]A", "YT[C>G]T", "YT[C>T]A", "YT[C>T]T"]
for qf in quin_files:

	# store sample ID from sample name
	id = ((qf[::-1].split("/")[0])[::-1]).replace("_SBS96_quintuplet_context.tsv", "")

	# read file and capture NT[C>N]NN mutations only
	qf_row = [id]
	qf_indiv = [(x.split("\t")[4][0])+(x.split("\t")[-1]) for x in (open(qf).read().rstrip("\n").split("\n")) if "T[C" in x.split("\t")[-1]]

	# add ambiguity R/Y
	ambig_qf_indiv = []
	for site in qf_indiv:
		if site.startswith(("A", "G")):
			ambig_qf_indiv += ["R" + site[1:]]
		elif site.startswith(("C", "T")):
			ambig_qf_indiv += ["Y" + site[1:]]

	# check for each pentamer
	for pt in apobec_pent_types:
		qf_row += [ambig_qf_indiv.count(pt)]

	qf_out += [qf_row]

# create dataframe
apobec_pent_counts = pd.DataFrame(qf_out, columns = (["ID"] + apobec_pent_types))
apobec_pent_counts.set_index("ID", inplace=True)
apobec_pent_counts.to_csv("APOBEC3A-vs-APOBEC3B_quintuplet_context.tsv", sep="\t", header=True, index=True)
