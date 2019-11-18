# Running PICRUSt on the UCT HPC cluster.

## This is based on using the ASV table produced by the dada2 pipeline here https://github.com/kviljoen/16S-rDNA-dada2-pipeline
1. Download ASV table created (e.g. 'ASV_counts.RDS'), import into R, standardize reads and merge at lowest available taxonomic level using the script 'ASV_de_novo_to_GG_mapping.R' (save phyloseq object as .RDS file for use in 3.)
#NOTE: You cannot run this script as is - you will need to adapt to specify your own files and filepaths.
#NOTE: You should run this script line-by-line not all in one go!
2. Install usearch11 (instructions https://www.drive5.com/usearch/manual/install.html). Run closed_ref_from_de_novo_forPICRUSt.sh with the fasta file with representative ASV seqs from step 1.
#NB: You cannot run this script as is - you will need to adapt to specify your own files and filepaths and set your HPC account group membership
3. You will need the .RDS file from 1. for the R script WISH_de_novo_to_closed_ref_OTUs.R together with the .txt mapping file created in 2. to create a ASV table with GG IDs and it's corresponding .biom file for use with PICRUSt.
4. Run PICRUSt.sh on .biom generated in 3 as follows:
4a. Log in to HPC
4b. Start an interactive job with something like
```
srun -N 1 --time=24:00:00 --pty bash
```
4c. Install picrust from bioconda: http://picrust.github.io/picrust/install.html#install
```
conda create -n picrust1 -c bioconda -c conda-forge picrust
conda activate picrust1
```
You may get an error that says your conda environment has not been initialized properly. If this happens do:

```
conda init bash
```
You should see different versions of python being used once you've activated the picrust1 environment.
```
(base) bash-4.2$ which python
/opt/exp_soft/anaconda3/bin/python
(base) bash-4.2$ conda activate picrust1
(picrust1) bash-4.2$ which python
/home/kviljoen/.conda/envs/picrust1/bin/python
```
You may get an error related to the python library numpy. If this happens you will have to install the h5py library (which has numpy as dependency) within your conda environment with:

```
conda install -c anaconda h5py
```
If you get an error stating disk quota exceeded please contact the system administrator to increase your chunk files.

4d. From the interactive job and activated picrust1 environment submit the PICRUSt.sh script BUT:
#NB: You cannot run this script as is - you will need to adapt to specify your own files and filepaths and set your HPC account group membership
```
sbatch PICRUSt.sh
```
