#!/bin/bash

#SBATCH -D /set/your/own/workdir/
#SBATCH -J de_novo_to_closed_ref
#SBATCH -o /set/your/own/workdir/.command.log
#SBATCH --no-requeue
#SBATCH -c 2
#SBATCH -t 05:00:00
#SBATCH --mem 4096
#SBATCH -p ada
#SBATCH --account cbio

inDir=/set/your/own/indir #specify input directory
outDir=/set/your/own/outdir #specify output directory

#specify your own copy of usearch11, input fasta and greengenes fasta file
/home/kviljoen/usearch11 -usearch_global $inDir/Levin_microbiome_merged_seqs.fasta -db /bb/DB/bio/qiime/greengenes/gg_13_8_otus/rep_set/97_otus.fasta  -id 0.97 -strand plus -uc $outDir/Levin_de_novo_repset_to_GG_13_8_map.txt 
#NOTE1: now download the resulting .txt file which maps de novo IDs to GG IDs and OTU rownames with GG IDs (for those that map) in R
