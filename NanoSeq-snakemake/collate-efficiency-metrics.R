library(dplyr)
library(stringr)
library(purrr)

# Efficiency estimation output files (one per sample)
files = list.files("12_efficiency_estimate", "*.tsv", full.names = T)

# Inititate list object
data_list <- list()

# Loop through files, extract contents, add as an element to list object, with sample name as the element name
for (i in 1:length(files)) { 
	data = read.table(files[i])
	colnames(data) <- c("stats", (str_split(files[i], "/", simplify = T)[,2] %>% str_split("\\.", simplify = T))[,1])
	data_list[[i]] <- data
}

# Convert list object to dataframe
efficiency <- reduce(data_list, full_join)

# Separate whole genome and chr1 read bundle metrics
wg_metrics <- efficiency[1:3,]
rb_metrics <- efficiency[4:11,]

# Save
write.table(rb_metrics, "12_efficiency_estimate/read-bundle-metrics-chr1-only.tsv", quote = F, row.names = F, sep = "\t")
write.table(wg_metrics, "12_efficiency_estimate/efficiency-metrics-whole-genome.tsv", quote = F, row.names = F, sep = "\t")

