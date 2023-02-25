#!/usr/bin/env bash

# Run from within data directory to get sbs mutation tables for nanoseq data 
# Usage: dbs_for_nanoseq_input_generator.sh  - run from data dir created using nanoflow workflow

# Note : idk if we need much of this at all, other than to gather the tables for output?

#If error, stop script
set -e

# Make output directory and change to it
mkdir sbs_files
cd sbs_files

# Copy sbs vcfs to new dir
echo "Copying SBS mutation files"
cp ../summary_PD*/*muts.tsv ./
echo "SBS mutation files copied!"

# I think we actually just need the above and then to output this directory to use in the R script?

# # Load farm version of bcftools
# module load bcftools

# # Filter vcfs to leave PASS variants only, and create new files with filtered only
# for file in *.vcf;
# do
# 	bcftools view -f PASS "$file" > temp_PASS.vcf
# 	mv temp_PASS.vcf "$(basename "$file" .vcf)"_PASS.vcf
# 	echo "$(basename "$file") filtered!"
# done

# # Delete extraneous files
# rm *.tbi
# rm *indel.vcf

# # Format vcf files into tables for DBS calling

# #Remove hashed lines at start
# for file in *.vcf;
# do
# 	grep -v "^#" "$file" > "$(basename "$file")"_table.tsv
# 	echo "$(basename "$file") converted to table!"
# done