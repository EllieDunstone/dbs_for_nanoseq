#!/usr/bin/env bash

# Input a canapps project number to get the zipped sbs vcfs for nanoseq data, unzip, delete tbi files, and filter for PASS only
# Then call dinucleotide mutations (using script based on caller from Inigo Martincorena, 2022)
# Usage: dinuc_nanoseq_caller_wrapper.sh  - run from data dir 

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

#Load R
module load R

Rscript /lustre/scratch119/casm/team154pc/ed4/phd/misc_scripts/dbs_calling_annotating_runs_subs.R



#echo "DBS called for $(basename "$file")!"






