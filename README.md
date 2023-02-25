# DBS for nanoseq

NB// WIP - no promises anything here works right now.

## Introduction

This repo contains scripts for the analysis of doublet base substitutions (DBS) using data generated using the NanoSeq pipeline at the Wellcome Sanger Institute. The scripts can be used to call DBS and other multi-nucleotide variants from the SBS calls output by the pipeline, generate mutation matrices for signature analysis, plot mutational spectra, and correct the raw mutation counts according to the dinucleotide coverage per sample.

## Scripts

The script "dbs_calling_annotating_runs_subs.R" annotates dinucleotides and runs of consecutive subs of any length using a function written by Iñigo Martincorena.

The script "dnv_caller.Rmd" generates a DBS78 matrix (see <https://cancer.sanger.ac.uk/signatures/dbs/>) for signature analysis from a table of DNVs. This is based on a script written by Andrew Lawson. You can generate the table of DNVs using the script "dbs_calling_annotating_runs_subs.R".

The script "dbs_call_correction.Rmd" uses the trinucleotide frequencies present in your NanoSeq data (output by the NanoSeq pipeline) to calculate the dinucleotide frequencies present, and then uses these values to correct your matrix of observed DBS counts according to the dinucleotide frequency in your data. This effectively maps your observed calls onto the frequency that would be expected in a normal diploid human genome.

The script "dbs_for_nanoseq_input_generator.sh" can be used by users of the NanoSeq pipeline at the Wellcome Sanger Institute to gather the mutation and ratio2genome files needed from the compute farm. This script assumes prior use of the nanoflow workflow to process the output data (not yet on github, will update).

Example input data is provided in the directory 'example_data'.


## Dependencies

The scripts load the following R libraries, which will need to be installed prior to running:

* tidyverse
* grid
* gtable
* deconstructSigs

## Usage

The script dbs_calling_annotating_runs_subs.R annotates dinucleotides and runs of consecutive subs of any length using a function written by Iñigo Martincorena. This function takes a data frame of mutations with the following columns: 
sampleID, chr = chromosome, pos = mutated site, ref = reference base, mut = mutant base.

etc etc add more here


## Acknowledgments

 Iñigo Martincorena and Andrew Lawson are thanked for their contributions to the code used in these scripts. Thank you to Federico Abascal for his help with implementing the correction of DBS counts according to dinucleotide frequency in the observed data, and writing of the original SBS code for NanoSeq that some of this code is based on.
