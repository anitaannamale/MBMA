# MBMA

A Mapping Based Microbiome Analysis tool for species and resistance genes quantification

## About
  A main challenge in the field of metagenomics is to quickly and precisely determine the composition and the relative abundance of constituents of a microbial community. State-of-the-art metagenomic methods rely on mapping-, kmer- or assembly- based  approaches for identification and quantification. However, currently available methods are either time-consuming, less sensitive and generates false positives, which leads to a lack of accuracy in quantification at gene and strain levels. 

  To overcome these limitations, we have implemented a tool called MBMA for Mapping Based Microbiome Analysis, a mapping-based approach for identification and quantification of species, strains and resistance genes, with three main innovations, the use of : <br />
- a **efficient and discriminatory database** for rapid quantification, <br />
- an **advanced counting method** to decrease the false discovery rate, <br />
- combined with **variant calling** from samples sequences for an accurate abundance prediction.

**MBMA** identifies and quantifies constituents (species, resistance genes) from metagenomic samples by mapping reads against a database and performing variant calling. Other constituents can be quantifies by changing the database. It has 3 way of working :
- **mapping**, it simply map reads against an indexed reference database, using different mapping tools (bowtie2, bwa and novoalign) and different counting methods (best, ex-aequo and shared), for non redundant databases. 
- **variant**, it performs, in addition to the mapping step, a variant calling step. Reads are mapped against a clustered database, and then variant calling using GATK is performed, to able an accurate quantification at a gene level.
- **mode**, it provide two presets modes for the identification of species (option : --species) and resistance genes (option : --resistance). For bacterial species, reads are mapped reads against *RefMG.v1*, and for resistance genes, against *ResFinder*. 

MBMA works with both Paired-end and Single-end reads. It can directly read input FASTQ files that are compressed using gzip. User should provide a database that is indexed with bowtie2, bwa or Novocraft.

MBMA run on executing tasks on a cluster. Informations about user's computing cluster should be given in the file *submission_template.sh*

## Version
1.0 

## Requirements

**Python** 2.7.8

One of those mapping tools : 

 **Bowtie2** 2.2.6 (require C++)
 **BWA** 0.7.7 (require C)
 **Novocraft** V3.02.12 (require C++)

For counting mapped reads : **pysam** 0.8.4

For variant calling :
**GenomeAnalysisTK** 3.4-0
**picard-tools** 1.94

**SAMtools** 1.2

## Download & Execute

To download MBMA clone the repository and execute, use the following commands :

```
git clone https://github.com/anitaannamale/MBMA.git
cd MBMA-master
./mbma 
```

## Usage



## Examples


## Code
This code is written in Python

File | Description
---|---
**count_matrix.py** | Code for merging several count tables into a counting matrix
**counting.py** | Creates counting tables by counting maped reads according to methods "best", "ex-aequo" or "shared"
**db_index_sub.sh** | Code for indexing the database for variant calling * *should be modified according to your computing cluster*
**help.py** | commandline parsing arguments
**mbma.py** | Main program
**submission_template.sh** | Code that launches task on the computing cluster * *should be modified according to your computing cluster*
**variant_predict.py** | Code for creating count table from the variant calling 

## Future work

## Test Datasets

## References
