source /local/gensoft2/adm/etc/profile.d/modules.sh # module location path

module add picard-tools/1.94 samtools/1.2 # add your modules 

# Create picard index if it doesn't exist
if [ ! -f {fasta_dict} ];
then
    echo "### Picard index creation "
    CreateSequenceDictionary R= {ref_fasta} O= {fasta_dict}
fi

# Create samtools index if it doesn't exist
if [ ! -f {fasta_fai} ];
then
    echo "### Samtools index creation "
    samtools faidx {ref_fasta}
fi
