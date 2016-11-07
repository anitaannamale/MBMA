
# MBMA

A Mapping Based Microbiome Analysis tool for species and resistance genes quantification

## About
  A main challenge in the field of metagenomics is to quickly and precisely determine the composition and the relative abundance of constituents of a microbial community. State-of-the-art metagenomic methods rely on mapping-, kmer- or assembly- based  approaches for identification and quantification. However, currently available methods are either time-consuming, less sensitive and generates false positives, which leads to a lack of accuracy in quantification at gene and strain levels. 

  To overcome these limitations, we have implemented a tool called MBMA for Mapping Based Microbiome Analysis, a mapping-based approach for identification and quantification of species, strains and resistance genes, with three main innovations, the use of : <br />
- a **efficient and discriminatory database** for rapid quantification, <br />
- an **advanced counting method** to decrease the false discovery rate, <br />
- combined with **variant calling** from samples sequences for an accurate abundance prediction.

**MBMA** identifies and quantifies constituents (species, resistance genes) from metagenomic samples by mapping reads against a database and performing variant calling. Other constituents can be quantifies by changing the database. It has 3 ways of working :
- **mapping**, it simply map reads against an indexed reference database, using different mapping tools (bowtie2, bwa and novoalign) and different counting methods (best, ex-aequo and shared), for non redundant databases. 
- **variant**, it performs, in addition to the mapping step, a variant calling step. Reads are mapped against a clustered database, and then variant calling using GATK is performed, to able an accurate quantification at a gene level.
- **mode**, it provide two presets modes for the identification of species (option : --species) and resistance genes (option : --resistance). For bacterial species, reads are mapped reads against *RefMG.v1*, and for resistance genes, against *ResFinder*. 

**MBMA** works with both Paired-end and Single-end reads. It can directly read input FASTQ files that are compressed using gzip. User should provide a database that is indexed with bowtie2, bwa or Novocraft.

**MBMA** works with nucleotide database 

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

To launch MBMA, computer cluster parameters should be specified in 3 files :
  - db_index_sub.sh : line 1 & 3
  - submission_template.sh  : line 1-14
  - job_make_matrices.sh : line 1-12

Modify the corresponding lines according to your computer cluster. 

**/!\ Do not modify the words in brackets**

## Usage
MBMA works with paired-end and single-end reads.

MBMA has 3 ways of working :
  - **mapping** : to map reads against a database and create a counting table. To launch MBMA mapping mode do `./MBMA mapping ` 
  - **variant** : to map reads against a database and perform variant calling using GATK. To launch MBMA variant mode do `./MBMA variant `
  - **mode** : is a preset of tools and database for mapping and variant modes. To launch MBMA mode do `./MBMA mode `

For more informations, see MBMA help, by running:
`/MBMA -h`, `MBMA mapping -h`, `MBMA variant -h`, or `MBMA mode -h` 

## MBMA mapping

Three mapping tool are available to map reads : **BWA**, **Bowtie2**, and **Novocraft**

**MBMA mapping** mode uses an **indexed database** with the selected mapping tool. 

Three counting methods are available to count read that map to a reference sequence : **best**, **ex-aequo**, and **shared**

Genes abundance is generated by two step. First, the unique mapped reads (reads mapping to a unique gene) were attributed to the corresponding genes. Then, the multiple reads (reads mapping to several genes) were attributed taking count of their weights depending on counting mode: best, ex-aequo and shared.  Abundance of a gene G (Ag)depends on the abundance of uniquely mapped reads (Au) and abundance of multiple reads (Am)  

![equation](https://latex.codecogs.com/gif.latex?A_%7Bg%7D%3D%20A_%7Bu%7D%20&plus;%20A_%7Bm%7D)

- **best** = a multiple read is attributed to a single reference gene randomly
- **ex-aequo** = a multiple read is attributed to each reference gene where it maps
- **shared** = a multiple read is attributed to each reference gene according to the ratio of their unique mapped counts called coefficient Co, as follow :

       ![equation](https://latex.codecogs.com/gif.latex?Co_%7Bi%7D%3D%20%5Cfrac%7BAu%7D%7B%5Csum_%7Bj%3D1%7D%5E%7BM%7D%20Au_%7Bj%7D%7D)
with i a multiple read, which maps on M genes.

       Thus, for a gene, the abundance of his multiple reads is :

       ![equation](https://latex.codecogs.com/gif.latex?A_%7Bm%7D%3D%20%7B%5Csum_%7Bi%3D1%7D%5E%7BN%7D%20Co_%7Bi%7D%7D)


    For example, if a gene G is mapped by 10 reads that only map to it (unique reads), and by one read that map to it and also on a gene X which is mapped by 5 unique reads, then :      ![equation](https://latex.codecogs.com/gif.latex?A_%7Bg%7D%3D%2010%20&plus;%20%5Cfrac%7B10%7D%7B10&plus;5%7D%5Capprox%2010.7)


**MBMA mapping** supports gzip compressed FASTQ input files and auto-detected from the file name.
You should provide a directory containing all the FASTQ files.
When using Paired-end data, paired files must be named : 

    file_R1.[fastq/fq] & file_R2.[fastq/fq]
    OR
    file_1.[fastq/fq] & file_2.[fastq/fq]
    

**MBMA ** provides the possibility to increase the mapping stringency level by using the option *--strict*. To increase the strengency of mapping, the minimum score threshold required for an alignement to be considered as "valid" is increased. For bowtie2, the strict score is 40 + 8.0 * ln(L), where L is the read length (default is 20 + 8.0 * ln(L)). For BWA, the strict score is 60 (default is 30). 

### Examples

    MBMA mapping -i data -o output -db database_index -e email@pasteur.fr -q queue --bowtie2 --shared
    MBMA mapping -i data -o output -db database_index -e email@pasteur.fr -q queue --bwa --best --strict
    MBMA mapping -i data -o output -db database_index -e email@pasteur.fr -q queue --bwa --ex-aequo --strict

## MBMA variant

**MBMA variant** is useful when the database is redundant, in cases where genes are too similar and may differ by only one nucleotide (for example resistant genes databases like ResFinder).
It allows an accurate quantification of genes for a redundant database.

**MBMA variant** works in two step process. First, reads are mapped against an reduced database, then using the alignement file (BAM file), variant calling is performed using GATK creating a variant profile, finally by comparing the variant profile with matrix of variant, the allele of genes are identified and their abundance are estimated.

A **reduced database of ResFinder** is given as test datasets. This reduced database is a clustered one by **CDHIT-est** at 98% of nucleotide identity in atleast 90% of coverage of the smallest sequence. From the formed clusters, a matrix of variant profile of the database's sequences is created. 

*How to create a own reduced database :*
User can create their own reduced database by running the script *make_matrices.py*

### Examples

    MBMMA variant -i data -o output -db database_index -fa database.fa -matrice VCF_matrices/ -t 4 -e email@pasteur.fr -q queue --bowtie2 --shared

## MBMA mode

**MBMA mode** is a preset set of options and arguments to identify and quantify bacterial species and resistant genes from metagenomics samples.

To run the identification and quantification of bacterial species from clinical samples, simply run the option **--species**:

    MBMA mode -i data -o output -t 8 -e email@pasteur.fr -q queue --species

This option will run MBMA mapping by aligning reads against RefMG.v1.padded database, using *bowtie2*, option *--strict* and *shared* counting method for bacterial species identification and quantification.

To run the identification and quantification of resistance genes from clinical samples, simply run the option **--resistance**:

    MBMA mode -i data -o output -e email@pasteur.fr -q queue --resistance

This option will run MBMA variant bu aligning reads against a reduced ResFinder database, using *bowtie2*, option *--strict*, *shared* counting method and provided *variant matrices* of ResFinder, for resistance genes identification and quantification. 

## Code
This code is written in Python. Are preceded by * executable code. Submission template code that should be modified according to your computing cluster are in *italic*

File | Description
---|---
**count_matrix.py** | Code for merging several count tables into a counting matrix
**counting.py** | Creates counting tables by counting maped reads according to methods "best", "ex-aequo" or "shared"
*db_index_sub.sh* | Submission script for indexing the database for variant calling * *should be modified according to your computing cluster*
**help.py** | commandline parsing arguments
* **mbma.py** | Main program
*submission_template.sh* | Code that launches task on the computing cluster * *should be modified according to your computing cluster*
**variant_predict.py** | Code for creating count table from the variant calling 
**cluster2multifasta.py** | Create a nucleotide sequence multifasta file for each cluster
*job_make_matrices.sh* | Submission script to clusterize a database and creates variant matrices. * *should be modified according to your computing cluster*
* **make_matrices.py** | Code that clusterize a database and creates variant matrices.

## Future work

## Test Datasets

Test datasets will be soon provided.

To test MBMA mapping (or MBMA mode --species) for bacterial species identification RefMG.v1.padded database is provided in the folder databases.

To test MBMA variant (or MBMA mode --resistance) for resistance gene identification the clustered ResFinder database is provided in the folder databases and his variant profile matrice is in the folder matrices.

## References

* Ciccarelli FD, Doerks T, von Mering C, Creevey CJ, Snel B, Bork P. (2006)
**Toward automatic reconstruction of a highly resolved tree of life.**
*Science* 311.5765 : 1283-1287. 
DOI:[10.1126/science.1123061](http://dx.doi.org/10.1126/science.1123061)

* Zankari E, Hasman H, Cosentino S, Vestergaard M, Rasmussen S, Lund O, Aarestrup FM, Larsen MV. (2012)
**Identification of acquired antimicrobial resistance genes.**
*Journal of antimicrobial chemotherapy*, 67.11: 2640-2644. 
DOI:[10.1093/jac/dks261](http://dx.doi.org/10.1093/jac/dks261)
