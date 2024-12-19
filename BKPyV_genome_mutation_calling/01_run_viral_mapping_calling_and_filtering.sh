#!/bin/sh
# ./01_run_viral_mapping_calling_and_filtering.sh genome.fa read1.fq.gz read2.fq.gz
# genome.fa should be faidx indexed (samtools)



## SCRIPT PURPOSE
# This script was used to call mutations of the BKPyV genome during infection of urothelial cells
# Reads are first mapped against the provided viral genome, corrected and mutations called
# Then custom filtering scripts are used to extract the SBS96 mutational signature catalogue (triplet and quintuplet
##



## ENVIRONMENT SETUP ##

# load modules
module load SAMtools 	## v1.20 used
module load minimap2	## v2.26 used
module load picard		## v2.20.2 used
module load python		## v3.11.3 used

# pip install pandas pysam collections --user 	# if needed

# define java settings (needed for many schedulers like slurm)
export MALLOC_ARENA_MAX=4
vmArgs="-Xmx36000M -XX:ParallelGCThreads=12"

# generate outname variable
outname=`echo $2 | rev | cut -d'/' -f1 | cut -d'_' -f2 | rev`



## PART 1 - map to viral genome ##

# use minimap2 to map to reference genome with -x sr for short read preset
minimap2 -ax sr -t 12 -2 $1 $2 $3 > $outname".sam"

# convert to BAM, keeping only mapped reads, sort and indexed
# use 12 threads with 2900MB RAM per thread (allow 3GB RAM per thread as samtools sort over-requests)
samtools view -b -F 4 -t $1.fai $outname".sam" > $outname".bam"
samtools sort -@ 12 -m 2900M -T $outname -o $outname".sorted.bam" $outname".bam"
samtools index $outname".sorted.bam"

# tidy up
rm $outname".sam" outname".bam"



## PART 2 - correct BAM and call variants ##

# use picard tools and samtools to correct the BAM file
# add read groups and sort by name for fixmate
java $vmArgs -jar picard.jar AddOrReplaceReadGroups I=$outname".sorted.bam" O=$outname".RG.bam" RGLB=$outname RGPU=$outname RGPL=illuminaNanoSeq RGSM=$outname
samtools index $outname".RG.bam"
samtools sort -n -@ 12 -m 2900M -T $outname -o $outname".RG.namesorted.bam" $outname".RG.bam"

# run fixmate
# -r removes secondary and unmapped
# -m adds mate score flags needed by markdup
samtools fixmate -r -m -@ 12 $outname".RG.namesorted.bam" $outname".RG.FM.bam"
samtools sort -@ 12 -m 2900M -T $outname -o $outname".RG.FM.sorted.bam" $outname".RG.FM.bam"
samtools index $outname".RG.FM.sorted.bam"

# tidy up
rm $outname".RG.namesorted.bam" $outname".RG.FM.bam" $outname".RG.bam" $outname".RG.bam.bai"

# run markdup
# -l is the read length
samtools markdup -l 150 -r -T $outname -S -@ 12 $outname".RG.FM.sorted.bam" $outname".RG.FM.MD.sorted.bam"
samtools index $outname".RG.FM.MD.sorted.bam"

# tidy up
rm $outname".RG.FM.sorted.bam" $outname".RG.FM.sorted.bam.bai"



## PART 3 - call variants using mpileup ##

# mpileup keeps base qual >=30 and mapping qual>=50, uses the ref genome to create a ref allele
# with mpileup indels are reported in a way that makes regex parsing difficult.
# A 3bp AGC insertion in a string of Ts would be reported: TTTTT+3AGCTTTT
# A 3bp TTT deletion from an AGTTTTAGTTTT repeating pattern would be reportted: AGTTTTAGT-3TTTAGTTTT
# The --no-output-ins and --no-output-del flags handle these (strangely). Using the flag once removes the indel sequence, but using twice also removes the indel size and symbol.
# mpileup -d 0 flag removes the max depth option
samtools mpileup -q 30 -Q 30 --no-output-ins --no-output-ins --no-output-del --no-output-del -f $1 -d 0 $outname".RG.FM.MD.sorted.bam" > $outname".MD.noindels.mpileup"


# need to process the mpileup files to get rid of reference bases, non ATGC characters etc etc
# capture cols 1-5 (don't care about qual etc as this has been filtered already)
# replace . and , and starting character and brakcets and $ use of asterisk - check the mpileup docs https://www.htslib.org/doc/samtools-mpileup.html
# make sure everything is uppercase
# keep lines where there is still some non-ref allele
# print with ratio length to check it isn't retained "germline" SNP

awk '$4>10 {print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5}' $outname".MD.noindels.mpileup" |
      sed 's/\.//g; s/\,//g; s/\^[a-zA-Z0-9]\?//g; s/]//g; s/\$//g; s/\*//g; s/@//g; s/\?//g; s/\[//g; s/\\//g' |
      tr '[:lower:]' '[:upper:]' |
      awk '$5!=""' |
      awk '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" (length($5)/$4) "\t" $5}' > $outname".accepted-variants.tsv"



## PART 4 - fitering and processing reading for SBS analysis

# process the accepted-variants into a form ready to assemble into an SBS catalogue
python 02_process_filtered_mpileup_variants.py $outname".accepted-variants.tsv"

# create the SBS96 catalogue file N[Ref>Alt]N
python 03_create_SBS96_catalogue.py $outname".accepted-variants.chr-pos-ref-alt.tsv" $1

# extract SBS96 mutations but with quintuplet context NN[Ref>Alt]NN
python 04_extract_SBS96_quintuplet_context.py $outname".accepted-variants.chr-pos-ref-alt.tsv" $1

## the APOBEC3A vs APOBEC3B pentamer dataframe is created across a whole set of quintuplet samples
## python 05_APOBEC3A-APOBEC3B_quintuplet_summary.py ./