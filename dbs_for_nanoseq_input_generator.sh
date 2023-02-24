#!/usr/bin/env bash

# Input a canapps project number to get the zipped sbs vcfs for nanoseq data, unzip, delete tbi files, and filter for PASS only
# Usage: dbs_for_nanoseq_input_generator.sh  - run from data dir 

#If error, stop script
set -e

# Make output directory and change to it
mkdir dinucs
cd dinucs

# Copy indel vcfs from irods
echo "Copying SBS mutation files"
#cp "/nfs/irods-cgp-*/intproj/${project_number}/sample/PD*/*muts.vcf*" ./
cp ../summary_PD*/*muts.tsv ./
echo "SBS mutation files copied!"

# Filter vcfs to leave PASS variants only, and create new files with filtered only
for file in *.vcf;
do
	bcftools view -f PASS "$file" > temp_PASS.vcf
	mv temp_PASS.vcf "$(basename "$file" .vcf)"_PASS.vcf
	echo "$(basename "$file") filtered!"
done

# Delete extraneous files
rm *.tbi
rm *indel.vcf

# Format vcf files into tables for DBS calling

#Remove hashed lines at start
for file in *.vcf;
do
	grep -v "^#" "$file" > "$(basename "$file")"_table.tsv
	echo "$(basename "$file") converted to table!"
done