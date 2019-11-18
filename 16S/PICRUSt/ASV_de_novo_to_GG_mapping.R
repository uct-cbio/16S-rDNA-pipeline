#SCRIPT PURPOSE: DEMONSTRATION OF 16S DOWNSTREAM ANALYSES (USING THE 16S PACKAGES metagenomeSeq, phyloseq, vegan)
#AUTHOR: KATIE LENNARD
#PLEASE SET WORKING DIRECTORY APPROPRIATELY

# Load custom functions, set working directories, import data --------------------------
source("microbiome_custom_functions.R")#Can be found here https://gist.github.com/kviljoen/97d36c689c5c9b9c39939c7a100720b9
#Read in ASV counts
ASVs <- readRDS("ASVs_counts.RDS"))#Default options with primer trimming & new refseq_rdp DB
str(ASVs)
dim(ASVs)
dim(ASVs)
ASVs <- t(ASVs)

#Read in taxonomy table
taxa <- readRDS("ASVs_taxonomy.RDS")
dim(taxa)

#Check ASV against taxonomy table for identical sample names in correct orientation
identical(rownames(ASVs),rownames(taxa))#TRUE

#Assign user-friendly ASV IDs to replace sequences e.g. ASV1, ASV2...
head(rownames(ASVs))
seqs <- rownames(ASVs)
ASV.IDs <- paste0("ASV",c(1:length(seqs)))
#Named vector:
names(seqs) <- ASV.IDs
head(seqs)

#Merge ASV table and taxonomic table as a phyloseq object
phy <- phyloseq(otu_table(ASVs,taxa_are_rows=TRUE),tax_table(taxa))
identical(taxa_names(phy),rownames(ASVs))#TRUE
taxa_names(phy) <- names(seqs)
str(phy)

#Let's explore the data
reads <- sample_sums(phy) #number of reads per sample
reads
#Let's standardize sample read count so we can compare between samples:
total = median(sample_sums(phy))#find median sample read count
standf = function(x, t=total) round(t * (x / sum(x)))#function to standardize to median sample read count
M.std = transform_sample_counts(phy, standf)#apply to phyloseq object
ntaxa(M.std)#number of taxa
sample_sums(M.std)#number of reads/sample

#----
#It's useful to collapse (sum) ASVs at their lowest available taxonomic annotations:
#For this we use the tax_glom.kv() function from the 'microbiome_custom_functions_tutorial.R' script loaded at the beginning of this script
M.phy <- tax_glom.kv(M.std)
## [1] "Removing phylogenetic tree"
## [1] "There are now xx merged taxa"
ntaxa(M.phy)
## [1] 122
head(tax_table(M.phy))
head(taxa_names(M.phy))
#1. We want a fasta file with ASV sequences to map against the greengenes database:
#But now - since we've merged ASVs - we have to select one representative ASV for each row/species. 
#Note that this selection is arbitrary
#You can check and see how much these seqs differ for a given species in a given sample. As an example:
P_copri <- grep("Prevotella_copri",as.character(unlist(tax_table(M.std)[,"Species"])))
tax_table(M.std)[P_copri,"Species"]
dim(tax_table(M.std)[P_copri,"Species"])
otu_table(M.std)[P_copri,]
copri_counts =data.frame(otu_table(M.std)[P_copri,])
dog16.test =which(copri_counts$Dog16!=0)
dog16_copri_ASVs = rownames(copri_counts)[dog16.test]
copri.16 = seqs[dog16_copri_ASVs]
library(ShortRead)#Install if you don't have this yet
writeFasta(ShortRead(sread = DNAStringSet(copri.16), id = BStringSet(names(copri.16))), 
           file = "dog_16_copri_ASVs_for_MAFFT.fasta")#You can use the online MAFFT platform to do multiple seq alignments
#----------------------
#So now lets get representative seqs for our merged ASVs
merged_seqs <- seqs[taxa_names(M.phy)]
head(merged_seqs)
length(merged_seqs)#
#Write to fasta file:
writeFasta(ShortRead(sread = DNAStringSet(merged_seqs), id = BStringSet(names(merged_seqs))), 
           file = "microbiome_merged_seqs.fasta")
#2. Also save ASV table with counts 
saveRDS(M.phy, file ="M_phy.RDS")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#PHASE 2
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#PURPOSE OF THIS SCRIPT: IMPORT THE RESULTS FROM MAPPING DE NOVO PICKED ASVs (USING THEIR SECQUENCES)
#AND USEARCH-GLOBAL TO GREENGENES IDS AND ASSIGN SUBSET THAT MATCHED TO GREENGENES WITH GREENGENES IDS (TYPICALLY THOSE WITH THE BEST
#GREENGENES ANNOTATIONS) TO THE ASV TABLE
#THE RESULTING ASV TABLE WITH GG IDS CAN NOW BE USED TO CONDUCT PICRUST ANALYSIS

#import mapping list of de novo IDs to GG IDs (This should have been produced on HPC with script closed_ref_from_de_novo_forPICRUSt.sh
o <- read.table(paste0(getwd(),"/de_novo_repset_to_GG_13_8_map.txt"), sep = "\t", header =F) 
head(o)#note: * means there was no GG ID match
#NB please check that the right columns are selected from 'o' for downstream use. When using the file de_novo_repset_to_GG_13_8_map.txt
#as is it should be columns 9 and 10 that are your de novo IDs and GG IDs respectively, but please check this.
length(which(o[,10]=="*"))#number of ASVs with no GG matches
o.closed <- o[o[,10]!="*",] #subset to exclude ASVs with no GG match
dim(o.closed)#number of OTUs that remain
head(o.closed)
rownames(o.closed) <- o.closed[,9]
#import ASV table from phase 1 (standardized, merged at lowest available taxonomic level)
M.phy <- load(paste0(getwd(),"M_phy.RDS"))#Change as appropriate - previously prepared ASV table from R, import as .RDS object 

#Extract ASV table
o.tab <- otu_table(M.phy)
dim(o.tab)#
length(which(rownames(o.tab)%in%o[,9]))#all should match?

#---------------------------------------------------
#filter asv table to exclude de novo asvs with no GG OTU matches
o.tab <- o.tab[rownames(o.tab)%in%o.closed[,9],]
dim(o.tab)#
#now substitute de novo IDs with GG IDs
o.closed <- o.closed[rownames(o.tab),]
rownames(o.tab) <- o.closed[,10]
head(o.tab)
dim(o.tab)#
#PROBLEM: how to deal with duplicate GG IDs - e.g. more than one of my de novo IDs map to the same GG ID
length(unique(rownames(o.tab)))#
length(which(duplicated(rownames(o.tab))))
#------
#example
rownames(o.tab)[85]
# rownames(o.tab)[85]
# [1] "198893"
# > cat("Synch1464614218233544000\n");
temp <- which(rownames(o.tab)=="198893")#84 85
o.tab[temp,]
x = o.closed[o.closed[,10]=="128382",]
x
# x
#         De_novo_ID GG_13_8_ID
# OTU_434    OTU_434     128382
# OTU_44      OTU_44     128382
# > cat("Synch1464614273695909000\n");
#end example
#------

#Will need to collapse ASVs with identical GG IDs somehow by summing the counts for those OTUs and replacing it with that total
#step 1. sort table by ID
o.t <- o.tab[ order(row.names(o.tab)), ]
head(rownames(o.t))#example output
# head(rownames(o.t))
# [1] "1000986" "1007750" "1028501" "1033552" "1035532" "1038074"
# > cat("Synch1464614289894893000\n");
#---
#step 2:
library(plyr)
t = data.frame(o.t)
t$ID <- rownames(o.t)
test=ddply(t, .(ID), function(x) colSums(x[,-dim(x)[2]], na.rm = TRUE))#break into chunks according to the column 'ID' and sum each chunk
head(test)
rownames(test) <- test$ID
dim(test)#
test <- test[,-1]
head(test)#tada!!!
#Let's convert to .biom format for use with PICRUSt
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomformat")
library(biomformat)
ASVs_biom <- make_biom(test, sample_metadata = NULL, observation_metadata = NULL,
          id = NULL, matrix_element_type = "int")
str(ASVs_biom)
#Now write biom file
write_biom(ASVs_biom, biom_file = "GG_IDs_forPICRUSt.biom")




