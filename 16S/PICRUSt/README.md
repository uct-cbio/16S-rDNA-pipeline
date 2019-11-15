# Running PICRUSt on the UCT HPC cluster.

##This is based on using the ASV table produced by the dada2 pipeline here https://github.com/kviljoen/16S-rDNA-dada2-pipeline
1. Download ASV table created (e.g. 'ASV_counts.RDS'), import into R, standardize reads and merge at lowest available taxonomic level using the script 'ASV_de_novo_to_GG_mapping.R' (save phyloseq object as .RDS file for use in 3.)
#NOTE: You cannot run this script as is - you will need to adapt to specify your own files and filepaths.
#NOTE: You should run this script line-by-line not all in one go!
2. Run closed_ref_from_de_novo_forPICRUSt.sh with the fasta file with representative ASV seqs from step 1.
3. You will need the .RDS file from 1. for the R script WISH_de_novo_to_closed_ref_OTUs.R together with the .txt mapping file created in 2. to create a ASV table with GG IDs and it's corresponding .biom file for use wit PICRUSt
4. Run PICRUSt.sh on output from 3.
