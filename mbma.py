#! /usr/bin/env python
# -*- coding: utf8 -*-

"""
mbma.py : MBMA - A Mapping Based Microbiome Analysis tool for
                species and resistance genes quantification
        
        This script, launch a qsub job to identify and quantify constituents 
        (species, resistance genes) from metagenomic samples by mapping reads 
        against a database and performing variant calling. It can be used to 
        quantify any other constituents by changing the database.
"""

__author__ = "Anita Annamalé"
__version__  = "1.0"
__copyright__ = "copyleft"
__date__ = "2016/07"

#-------------------------- MODULES IMPORTATION -------------------------------#

import sys
import os
import argparse
import re
import subprocess
import shlex
import glob
import time
import shutil
import stat

import help as h


#-------------------------- FUNCTIONS DEFINITION ------------------------------#

def is_database(dbname, prog):
    """
    Check if the database is indexed with the mapping tool
    
    Args:
        dbname [string] = database name
        prog [string] = expected extension of the database of the mapping tool

    Returns :
        database [string] = for NovoAlign, the path to the file,
                            for Bowtie2 and BWA, the path to the database 
                                directory 
        or quit
    """

    db = os.path.abspath(dbname) # absolute path of the database

    if prog == 'novo':
        if not os.path.isfile(db):
            sys.exit("[ERROR] The database isn't indexed with NovoAlign")
        return db

    db_dir = os.path.dirname(db) # database directory
    db_basename = os.path.splitext(os.path.basename(db))[0]

    # for bowtie2, check 6 files exists in the directory with the index 
    # extension .bt2
    if prog == 'bowtie2':
        ext = ['.bt2', '.bt2l']
        onlyfiles = [f for f in os.listdir(db_dir) 
                     if os.path.isfile(os.path.join(db_dir, f))
                     and f.startswith(db_basename) 
                     and (os.path.splitext(f)[1] in ext)]
        if not len(onlyfiles) == 6:
            sys.exit("[ERROR] The database isn't indexed with Bowtie2")
        return db

    # for bwa, check if files with the expected extension exists
    if prog == 'bwa':
        ext = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        onlyfiles = [f for f in os.listdir(db_dir) 
                     if os.path.isfile(os.path.join(db_dir, f))
                     and f.startswith(db_basename) 
                     and (os.path.splitext(f)[1] in ext) ]
        if not len(onlyfiles) == 5:
            sys.exit("[ERROR] The database isn't indexed with BWA")
        return db


def create_index(filename, text, queue, index_name):
    """

    """
    # create a submission file for qrsh 
    with open(filename, "wt") as fileout:
        fileout.write(text)

    # make an executable file
    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IEXEC)

    # launch the qrsh
    commandline = "qrsh -q {0} -cwd ./{1}".format(queue, filename)
    arg_split = shlex.split(commandline)
    index = os.path.splitext(index_name)[1]
    print "Variant calling : Creating {0} index".format(index[1:].upper())  
    try:
        p = subprocess.call(arg_split)
    except subprocess.CalledProcessError as e:
        p = e.p
        sys.exit(e.returncode)

    # wait until the index is created and then delete the file
    while not os.path.exists(index_name):
        time.sleep(1)
    os.remove(filename)


def find_fastq_files(file_list, mode):
    """
    Find fastq file by extension according to mode (Paired-end or Single end) 
    and rename them.
    
    Args:
        file_list [list] =  all files of a directory
        mode [string] = 'Single End' or 'Paire End'
    
    Returns:
        ext [str] = extension of FASTQ files
        num [int] = number of task to execute
        files_name [list] = list of all fastq filenames
    """

    ext_list = []
    filenames = []
    
    # Get a list of fastq files
    for file1 in file_list:
        split_filename = os.path.splitext(file1)
        ext_1 = split_filename[1]
        
        # check if compressed file
        if ext_1 == ".gz":
            ext = '.gz'
            split_filename = os.path.splitext(split_filename[0])
            ext_1 = split_filename[1]
        else:
            ext = ""
            
        # get only fastq files
        if (ext_1 == '.fq') or (ext_1 == '.fastq'):
            ext = ext_1 + ext
            if mode == 'SE':
                filenames.append(file1)
            # if mode PE
            else : 
                name_wo_ext = os.path.splitext(split_filename[0])[0]
                if (name_wo_ext[-1] == '1'):
                    mate_file = re.sub('1' + ext, '2' + ext, file1)
                    # check if mate is present
                    if mate_file in file_list:
                        filenames.append([file1, mate_file])
                    ext = "1" + ext
                    # get filename format
                    if (name_wo_ext[-2] == 'R') : 
                        ext = 'R' + ext
                    if (name_wo_ext[-2] == '_' or name_wo_ext[-3] == '_'):
                        ext = "_" + ext
                    ext_list.append(ext)

    # Get the number of task to execute
    num = len(filenames)
    if num == 0:
        sys.exit("[FATAL] No FASTQ file detected. FASTQ file names must have "
                 "the format: file_R1.fastq & file_R2.fastq. OR file_1.fastq "
                 "& file_2.fastq")

    # Get the extension of fastq files
    if mode == 'PE': 
        # check if all fastq files have the same extension
        if all(ext==ext_list[0] for ext in ext_list):
            ext = ext_list[0]
        else :
            sys.exit("[FATAL] All FASTQ files must have the same name format. "
                   "Example for Paired-end reads: "
                   "file_[R1/1].[fastq/fq] & file_[R2/2].[fastq/fq].")
    return ext, num, filenames


def create_dir(name, variant = False):
    """
    Function that create directories in the output directory: 
        - err, out, sam, bam, comptage for mapping mode
        - in addition: vcf and new_comptage for variant calling mode

    Args:
        name [string] = output directory name
        variant [True/False] = optinal argument, set to True for variant 
                               calling mode, [default : False]
    """
    try:
        os.mkdir(name)
        os.mkdir("{0}/out".format(name))
        os.mkdir("{0}/sam".format(name))
        os.mkdir("{0}/bam".format(name))
        os.mkdir("{0}/comptage".format(name))
        if variant: # if variant calling mode
            os.mkdir("{0}/vcf".format(name))
            os.mkdir("{0}/new_comptage".format(name))

    except OSError as ose:
        sys.exit("Cannot create {0} directory, "
                 "it already exists".format(ose.filename))


def write_tmp_file(file_list, output):
    """
    Write a file containing all the FASTQ file. This file will be read
    by the qsub submission script.

    Args:
        file_list [list] = list of fastq filenames
        output [string] = the output filename

    Returns:
        output [string] = output filename
    """
    output = "{0}/Files2read.txt".format(args['outdir']) # output name
    
    # write file
    with open(output, "wt") as out:
        for element in file_list:
            if len(element) == 2:
                out.write('\t'.join(element) + '\n')
            else :
                out.write(element + '\n')

    return output


def write_soum(template, param):
    """
    Write the qsub submission script.

    Args:
        param [dict] = containing all information needed for the submission 
                       script creation
    
    No returns
    """
    # read the template
    filein = open(template, "rt")
    temp = filein.read()
    # Write the submission script
    with open("{0}/submission.sh".format(param['outdir']), "wt") as out:
        out.write(temp.format(**param))
    # close the template
    filein.close()


def progress(count, total, suffix=''):
    """
    Function that prints a progress bar.

    Args:
        count [int] = number of executed events
        total [int] = total number of events
        suffix [str] = supplementary text to print
    """
    bar_len = 60
    filled_len = int(round(bar_len * count/float(total)))
    
    percents = round(100.0 * count / float(total),1)
    bar = '=' * filled_len + '-' * (bar_len - filled_len)

    sys.stdout.write('[%s] %s%s ..%s\r' % (bar, percents, '%', suffix))
    sys.stdout.flush()


#--------------------------------- MAIN ---------------------------------------#

if __name__ == '__main__' :
    
    # Creating mbma parser
    arg_parser = h.mbma_parser()

    # print help
    if not len(sys.argv) >= 2 :
        arg_parser.print_help()
        exit(1)

    # locate script dir
    script_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    
    ## PARSING COMMANDLINE -----------------------------------------------------

    # Parse commandline arguments
    try:
        args = arg_parser.parse_args()
        args = vars(args)
    except:
        exit(1)

    # check threads argument    
    if args['threads'] < 1:
        sys.stderr.write("[Warning] Illegal number of threads provided: %i, "
                          "MBMA will use a single thread by default" 
                          % args['threads'])
        args['threads'] = 1

    # cleaning the dictionnary
    for key, value in args.items():
        if value == None:
            del(args[key])

        # working modes parameters
    if 'species' in args:
        args.update({'bowtie2':'True',
                     'shared':'True',
                     'strict':'strict',
                     'ref_db':script_path + '/databases/index/refmg_bowtie2/' + 
                              'refmg'})

    if 'resistance' in args:
        args.update({'bowtie2':'True',
                     'shared':'True',
                     'strict':'strict',
                     'ref_db':script_path + '/databases/index/' + 
                              'clustered_ResFinder_bowtie2/resfinder98',
                     'ref_fasta':script_path + '/databases/' + 
                                 'clustered_ResFinder.fa',
                     'variant_matrice':script_path + '/matrices/' + 
                                       'clustered_ResFinder/'})



    # check database's index
        
        # mapping tool
    if 'bowtie2' in args:
        args['ref_db'] = is_database(args['ref_db'], 'bowtie2')
    elif 'bwa' in args:
        args['ref_db'] = is_database(args['ref_db'], 'bwa')
    elif 'novo' in args:
        args['ref_db'] = is_database(args['ref_db'], 'novo')
        
        # variant calling
    if 'ref_fasta' in args:
        fasta = args['ref_fasta']
        fasta_dict = os.path.splitext(fasta)[0] + ".dict"
        fasta_fai = fasta + '.fai'
        # picard index checking
        if not os.path.isfile(fasta_dict):
            filename = 'dict_index.sh'
            text = ("source /local/gensoft2/adm/etc/profile.d/modules.sh\n"
                    "module add picard-tools/1.94\n"
                    "CreateSequenceDictionary R= {0} O= {1}\n".format(fasta, 
                                                                    fasta_dict))
            create_index(filename, text, args['queue'],fasta_dict)

        # samtools index checking
        if not os.path.isfile(fasta_fai):
            filename = 'fai_index.sh'
            text = ("source /local/gensoft2/adm/etc/profile.d/modules.sh\n"
                    "module add samtools/1.2\n"
                    "samtools faidx {0}\n".format(fasta))
            create_index(filename, text, args['queue'],fasta_fai)


    ### INPUT ------------------------------------------------------------------
    
    files = [] # files to process
    
    # get absolute path of all files in input directory
    files = os.listdir(args['indir'])
    files = [ '{0}/{1}'.format(args['indir'],file1) for file1 in files]
    
    # number of tasks
    args['ext'], args['num'], files_name = find_fastq_files(files, args['mode'])


    ### OUTPUT -----------------------------------------------------------------

    # create output directories
    if 'ref_fasta' in args:
        create_dir(args['outdir'], variant=True)
    else:
        create_dir(args['outdir'])
    
    # create a temporary file
    args['read_file'] = write_tmp_file(files_name, args['outdir'])

    
    # WRITE QSUB SOUMISSION SCRIPT ---------------------------------------------

    # input files
    if (args['mode'] == 'PE'):
        args['input'] = ('INFILE_1=`sed "${{SGE_TASK_ID}}q;d" {read_file} |'
                         ' cut -f1`\n'
                         'INFILE_2=`sed "${{SGE_TASK_ID}}q;d" {read_file} | '
                         'cut -f2`'.format(**args))
    else :
        args['input'] = ('INFILE_1=`sed "${{SGE_TASK_ID}}q;d" '
                         '{read_file}`'.format(**args))

    # output directory
    args['stdout'] = '{0}/out'.format(args['outdir'])

    # mapping tools    
        # Bowtie2
    if 'bowtie2' in args :
        args['map_tool'] = args['bowtie2']
        
        # commandline for bowtie2
        cmd = ("bowtie2 -p {threads} -x {ref_db} -q --local "
               "--sensitive-local ".format(**args))

        # input
        if args['mode'] == 'PE':
            cmd += "-1 $INFILE_1 -2 $INFILE_2 "
        else :
            cmd += "-U ${INFILE_1} "
        
        # mapping strict
        if 'strict' in args:
            cmd += "--score-min G,40,8 "

        # get all alignements
        if 'ex_aequo' in args or 'shared' in args :
            cmd += "-a "

        # output
        cmd += "-S $OUTDIR/sam/${BASE}.sam "
        
        # get cmd 
        args['map_cmd'] = cmd
    
        ### BWA
    elif 'bwa' in args :
        args['map_tool'] = args['bwa']
        
        # commandline for bwa
        cmd = "bwa mem -t {threads} {ref_db} ".format(**args)
        
        # input
        if args['mode'] == 'PE':
            cmd += "$INFILE_1 $INFILE_2 "
        else :
            cmd += "$INFILE_1 "

        # mapping strict
        if 'strict' in args:
            cmd += "-T 60 "

        # get all alignements
        if 'ex_aequo' in args or 'shared' in args :
            cmd += "-a "

        # output
        cmd +=  "> $OUTDIR/sam/${BASE}.sam "

        # get cmd
        args['map_cmd'] = cmd
    
        # Novoalign
    elif 'novo' in args :
        args['map_tool'] = args['novo']
        
        # commandline for novo
        cmd = "novoalign -d {ref_db} -F STDFQ -o SAM ".format(**args)
        
        # input
        if args['mode'] == 'PE':
            cmd += "-f $INFILE_1 $INFILE_2  "

        # get all alignements
        if 'ex_aequo' in args or 'shared' in args :
            cmd += "-r All "
        else :
            cmd += '-r Random '

        # output
        cmd +=  "> $OUTDIR/sam/original_${BASE}.sam \n"

        # clean sam file, by removing bad alignements (score = 0)
        cmd += ("grep -v 'AS:i:0' $OUTDIR/sam/original_${BASE}.sam "
                "> $OUTDIR/sam/${BASE}.sam ")
        
        # get cmd
        args['map_cmd'] = cmd
    
    # mapping tools
    args['script_loc'] = script_path
    if 'best' in args:
        args['count_tool'] = 'best'
    elif 'ex_aequo' in args:
        args['count_tool'] = 'ex_aequo'
    elif 'shared' in args:
        args['count_tool'] = 'shared'
    if 'ref_fasta' in args:
        args['bam'] = '-bam'
    else :
        args['bam'] = ''
   
    # variant calling
    args['variant']= ""
    if 'ref_fasta' in args:
        args['variant'] = ("""
start_time_variant=$(timer)
echo "Variant calling analysis with GATK started at $(date +"%T")"
java -jar  \
/local/gensoft2/exe/picard-tools/1.78/libexec/AddOrReplaceReadGroups.jar \
I=$OUTDIR/bam/filtered_${{BASE}}.bam O=$OUTDIR/bam/grp_${{BASE}}.bam RGID=4 \
RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20

samtools index $OUTDIR/bam/grp_${{BASE}}.bam

GenomeAnalysisTK -T SplitNCigarReads -R {ref_fasta} \
-I $OUTDIR/bam/grp_${{BASE}}.bam -o $OUTDIR/bam/bon_${{BASE}}.bam \
-rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS

GenomeAnalysisTK -T HaplotypeCaller -R {ref_fasta} \
-I $OUTDIR/bam/bon_${{BASE}}.bam --genotyping_mode DISCOVERY \
-stand_call_conf 0 -stand_emit_conf 0 -o $OUTDIR/vcf/${{BASE}}.vcf

python {script_loc}/variant_predict.py {variant_matrice} \
$OUTDIR/vcf/${{BASE}}.vcf $OUTDIR/comptage/${{BASE}}.txt \
$OUTDIR/new_comptage/${{BASE}}.txt
echo "Elapsed time with GATK : $(timer $start_time_variant)
""".format(**args))


    # LAUNCH QSUB SOUMISSION SCRIPT --------------------------------------------

    submission_template = script_path + '/submission_template.sh'
    try :
        write_soum(submission_template, args)
    except (IOError, OSError) :
        shutil.rmtree(args['outdir'])
        sys.exit("Error : Writing submission script")
    except :
        shutil.rmtree(args['outdir'])
        sys.exit("Unexpected error: {0}".format(sys.exc_info()[0]))

    commandline = "qsub {outdir}/submission.sh".format(**args)
    arg_split = shlex.split(commandline)
    
    with open("{0}/stdout.txt".format(args['outdir']),"at") as out, \
         open("{0}/stderr.txt".format(args['outdir']),"at") as err:
        try:
            p = subprocess.call(arg_split, stdout=out, stderr=err)

        except subprocess.CalledProcessError as e:
            p = e.p
            sys.exit(e.returncode)

    errfile = "{outdir}/stderr.txt".format(**args)
    if os.stat(errfile).st_size != 0:
        with open(errfile, "at") as errfile:
            sys.stderr.write(errfile.read())
        
        shutil.rmtree(args['outdir'], ignore_errors=True)
        sys.exit(1)
    
    
    # CHECK PROGRESSION --------------------------------------------------------
    
    nb_tables = 0
        
    while (nb_tables < args['num']):
        time.sleep(2)
        # directory of the counting tables
        count_table_dir = 'comptage'
        if 'ref_fasta' in args:
            count_table_dir = 'new_comptage'
        # getnb of created files in the directory 
        nb_tables = len(glob.glob('{0}/{1}/*.txt'.format(args['outdir'], 
                                                         count_table_dir)))
        # print progression
        progress(nb_tables, args['num'], 
                '{0} count table out of {num} created'.format(nb_tables,**args))


    # COUNT MATRIX CREATION ----------------------------------------------------
    
    if (nb_tables == args['num']):
        
        # delete temporary file
        os.remove(args['read_file'])

        # delete old count table directory and rename the new as "comptage" :
        if 'ref_fasta' in args:
            old_dir = args['outdir'] + '/comptage'
            new_dir = args['outdir'] + '/new_comptage'
            shutil.rmtree(old_dir)
            os.rename(new_dir, old_dir)

        # count matrix creation
        with open("{0}/stdout.txt".format(args['outdir']),"at") as out, \
             open("{0}/stderr.txt".format(args['outdir']),"at") as err:
            try :
                commandline = ("python {0}/count_matrix.py -d {1}/comptage "
                               "-o {1}/comptage/count_matrix.txt".format(script_path,
                                                                args['outdir']))
            
                args = shlex.split(commandline)
                p = subprocess.call(args,stdout=out, stderr=err)
            
            except subprocess.CalledProcessError as e:
                p = e.p
                sys.exit(e.returncode)

    sys.exit("\n\nMBMA successfully completed")