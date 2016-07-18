# MBMA

A Mapping Based Microbiome Analysis tool for species and resistance genes quantification

## About
  A main challenge in the field of metagenomics is to quickly and precisely determine the composition and the relative abundance of constituents of a microbial community. State-of-the-art metagenomic methods rely on mapping-, kmer- or assembly- based  approaches for identification and quantification. However, currently available methods are either time-consuming, less sensitive and generates false positives, which leads to a lack of accuracy in quantification at gene and strain levels. 

  To overcome these limitations, we developed Mapping Based Microbiome Analysis (MBMA), a mapping-based approach for identification and quantification of species, strains and resistance genes, with three main innovations, the use of : <br />
- a **efficient and discriminatory database** for rapid quantification, <br />
- an **advanced counting method** to decrease the false discovery rate, <br />
- combined with **variant calling** from samples sequences for an accurate abundance prediction.

**MBMA** identifies and quantifies constituents (species, resistance genes) from metagenomic samples by mapping reads against a database and performing variant calling. Other constituents can be quantifies by changing the database. It has 3 way of working :
- **mapping**, it simply map reads against an indexed reference database, using different mapping tools (bowtie2, bwa and novoalign) and different counting methods (best, ex-aequo and shared), for non redundant databases. 
- **variant**, it performs, in addition to the mapping step, a variant calling step. Reads are mapped against a clustered database, and then variant calling using GATK is performed, to able an accurate quantification at a gene level.
- **mode**, it provide two presets modes for the identification of species (option : --species) and resistance genes (option : --resistance). For bacterial species, reads are mapped reads against *RefMG.v1*, and for resistance genes, against *ResFinder*. 

## Version
1.0

## Requirements

## Usage


## Examples


