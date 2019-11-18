#!/bin/bash

#SBATCH -D /specify/your/own/working/directory
#SBATCH -J picrust
#SBATCH -o /specify/errorlog/.picrust.log
#SBATCH --no-requeue
#SBATCH -c 4
#SBATCH -t 05:00:00
#SBATCH --mem 8192
#SBATCH -p ada
#SBATCH --account cbio

inDir=/specify/you/own/input/directory
outDir=$inDir #specify output directory

normalize_by_copy_number.py -i $inDir/GG_IDs_forPICRUSt.biom -o $outDir/CN_norm.biom

predict_metagenomes.py -i $outDir/CN_norm.biom -o $outDir/predicted_metagenome.biom --accuracy_metrics $outDir/predicted_metagenome.NSTIs.txt

#with_confidence returns metagenomic 95% confidence intervals for each gene prediction
#the --accuracy_metrics flags returns NSTI scores (The NSTI scores that you get back are the branch length on the Greengenes tree separating taxa to be predicted from the nearest sequenced genome. This set of these distances is averaged, weighting by the abundance of each organism)
#Note: used GG13.8 to map our IDs to greengenes but PICRUSt uses 13.5. This appears to be fine however as stated by one of the PICRUSt authors Daniel Mcdonald here https://groups.google.com/forum/#!topic/picrust-users/LTxmapB1xiA

categorize_by_function.py -i $outDir/predicted_metagenome.biom -l 3 -c KEGG_Pathways -o $outDir/categorize_by_function_level_3.biom
categorize_by_function.py -i $outDir/predicted_metagenome.biom -l 2 -c KEGG_Pathways -o $outDir/categorize_by_function_level_2.biom
categorize_by_function.py -i $outDir/predicted_metagenome.biom -l 1 -c KEGG_Pathways -o $outDir/categorize_by_function_level_1.biom

#convert .biom files to .txt files for use in R
biom convert -i $outDir/predicted_metagenome.biom -o $outDir/predicted_metagenome.txt -b #b specifies biom --> txt
biom convert -i $outDir/categorize_by_function_level_3.biom -o $outDir/categorize_by_function_level_3.txt -b
biom convert -i $outDir/categorize_by_function_level_2.biom -o $outDir/categorize_by_function_level_2.txt -b
biom convert -i $outDir/categorize_by_function_level_1.biom -o $outDir/categorize_by_function_level_1.txt -b

#STAMP software can also use .txt files the headers just need minor modifications I think (see exisiting .txt files previously used for STAMP) 

#lets see which OTUs underly the predicted KEGG functions using metagenome_contributions.py
metagenome_contributions.py -i $outDir/CN_norm.biom -o $outDir/metagenome_contributions.txt -t ko

