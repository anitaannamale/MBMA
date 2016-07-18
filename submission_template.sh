#!/bin/sh

#$ -S /bin/bash
#$ -M {email}
#$ -m bea
#$ -q {queue}
#$ -N {jobname}
#$ -j y
#$ -o {stdout}
#$ -pe thread {threads} 
#$ -t 1-{num}

### Import your modules
module add Python/2.7.8 bowtie2/2.2.6 bwa/0.7.7 novocraft/V3.02.12 samtools/1.2 GenomeAnalysisTK/3.4-0 picard-tools/1.94 

### INPUT FILES
{input}

### INPUT FILE BASENAME
BASE=`basename ${{INFILE_1%{ext}}}`

### OUTPUT DIRECTORY
OUTDIR={outdir}

### FUNCTIONS
function timer()
{{
    if [[ $# -eq 0 ]]; then
        echo $(date '+%s')
    else
        local  stime=$1
        etime=$(date '+%s')
        if [[ -z '$stime' ]]; then stime=$etime; fi
        dt=$((etime - stime))
        ds=$((dt % 60))
        dm=$(((dt / 60) % 60))
        dh=$((dt / 3600))
        printf '%d:%02d:%02d' $dh $dm $ds
    fi
}}

###PROGRAMS
export LD_LIBRARY_PATH=/pasteur/projets/Matrix/metagenomics/htslib/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/pasteur/projets/Matrix/metagenomics/python-lib/lib/python2.7/site-packages/:$PYTHONPATH

echo $date
echo "Mapping files: $INFILE_1 & $INFILE_2" 
start_time=$(timer)
echo "Started at " $(date +"%T") 

start_time_{map_tool}=$(timer)
echo "Mapping with {map_tool} started at $(date +"%T")"
{map_cmd}
echo "Elapsed time with {map_tool} : $(timer $start_time_{map_tool})"

start_time_{count_tool}=$(timer)
echo "Mapping analysis with {count_tool} started at $(date +"%T")"
python {script_loc}/counting.py $OUTDIR/sam/${{BASE}}.sam \
$OUTDIR/comptage/${{BASE}}.txt --{count_tool} {bam}
echo "Elapsed time with {count_tool} : $(timer $start_time_{count_tool})"
# if variant calling mode
{variant}
echo "Total duration :  $(timer $start_time)"