#!/bin/sh

#$ -S /bin/bash
#$ -M {email}
#$ -m bea
#$ -q {queue}
#$ -N {jobname}
#$ -j y
#$ -o {stdout}

### Import your modules
module add Python/2.7.8 java/1.8.0 cd-hit/v4.6.1 muscle/3.8.31 

###PROGRAMS
# clustering
cd-hit-est -i {database} -c {id_threshold} -aS 0.9 -d 0 -g 1 -r 1 -o {clustered_db}

# create multifasta
python {script_loc}/cluster2multifasta.py {clustered_db}.clstr {database} {clstr_fasta}

# filter only multifasta
mkdir {clstr_fasta}/cluster
for file in {clstr_fasta}/*.fa; do a=`grep '>' $file | wc -l`; if [ ! $a == 1 ] ; then mv $file {clstr_fasta}/cluster/ ; fi; done;

# alignement
mkdir {output}/alignment
for file in {clstr_fasta}/cluster/*.fa; do muscle -in $file -out {output}/alignment/`basename $file` -clwstrict; done;
for file in {output}/alignment/*.fa; do mv $file  ${{file%fa}}aln; done;

# VCF calling
mkdir {output}/vcf_db
for file in {output}/alignment/*.aln ; do {script_loc}/tools/jvarkit/dist/msa2vcf $file -c `basename ${{file%.aln}}` -R `basename ${{file%.aln}}` -o {output}/vcf_db/`basename ${{file%aln}}`vcf ; done;

