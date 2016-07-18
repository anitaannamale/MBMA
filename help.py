#! /usr/bin/env python
# -*- coding: utf8 -*-

"""
help.py : MBMA - A Mapping Based Microbiome Analysis tool for
                species and resistance genes quantification
        
        Module containing all functions to parse commandline arguments
"""

__author__ = "Anita Annamalé"
__version__  = "1.0"
__copyright__ = "copyleft"
__date__ = "2016/07"

#-------------------------- MODULES IMPORTATION -------------------------------#

import sys
import os
import argparse
from argparse import RawTextHelpFormatter

#---------------------------- CLASS DEFINITION --------------------------------#

class color:
   BOLD = '\033[1m'
   END = '\033[0m'
   
#-------------------------- FUNCTIONS DEFINITION ------------------------------#

def is_dir(path):
    """
    Check if a given path is an existing directory.

    Args:
        path [string] = Path to the directory

    Returns:
        abs_path [string] = absolute path
        or quit
    """
    abs_path = os.path.abspath(path) # get absolute path

    # if not a directory
    if not os.path.isdir(abs_path):
        # check if it is a path to a file
        if os.path.isfile(abs_path):
            msg = "{0} is a file not a directory.".format(abs_path)
        # else the path doesn't not exist
        else:
            msg = "The path {0} does not exist.".format(abs_path)
        raise argparse.ArgumentTypeError(msg)

    return abs_path


def dir_not_exists(path):
    """
    Check if a path exist or not.

    Args:
        path [string] = Path

    Returns:
        abs_path [string] = absolute path if it doesn't exist
        or quit
    """
    abs_path = os.path.abspath(path)

    # if path exists
    if os.path.exists(abs_path):
        msg = ("The path {0} already exist. "
               "Please give another directory name".format(abs_path))
        raise argparse.ArgumentTypeError(msg)

    return abs_path


def is_mode(text):
    """
    Check if the sequencing mode is 'PE' or 'SE', for reads.

    Args:
        text [string] = the text to check

    Returns:
        clean_text [string]
    """
    clean_text = text.upper()
    
    # if the text is not SE or PE
    if not (clean_text == 'SE' or clean_text == 'PE'):
        msg = "{0} isn't a valid mode ['PE' or 'SE' only].".format(clean_text)
        raise argparse.ArgumentTypeError(msg)

    return clean_text


def is_fasta(path):
    """
    Check if a path is an existing fasta file, only checks the extension of the 
    file.

    Args:
        path [string] = Path to the file

    Returns:
        abs_path [string] = absolute path to the fasta file if it exists
        or quit
    """
    abs_path = os.path.abspath(path) # get absolute path

    # if not a file
    if not os.path.isfile(abs_path):
        # check if it is a path to a directory
        if os.path.isdir(abs_path):
            msg = "{0} is a directory not a file.".format(abs_path)
        # else the path doesn't not exist
        else:
            msg = "The path {0} does not exist.".format(abs_path)
        raise argparse.ArgumentTypeError(msg)

    # else, it's a file, check the extension
    else:
        split_filename = os.path.splitext(abs_path)
        ext = split_filename[1] # get the extension
        if not (ext == '.fa') or (ext == ".fasta"):
            msg = ("The given file '{0}' isn't a fasta file, "
                   "please provide a fasta file [.fa/.fasta]".format(abs_path))
            raise argparse.ArgumentTypeError(msg)

    return abs_path


# HELP for MBMA

def mbma_parser():
    """
    MBMA help and arguments
    
    No arguments
    
    Return : parser
    """
    parser = argparse.ArgumentParser(
        formatter_class = RawTextHelpFormatter,
        
        description = color.BOLD + '\r       ' +
"""
             MBMA - A Mapping Based Microbiome Analysis tool for
                species and resistance genes quantification  

DESCRIPTION

    MBMA """ + color.END + """ identifies and quantifies constituents (species, resistance genes) 
    from metagenomic samples by mapping reads against a database and performing 
    variant calling. It can be used to quantify any other constituents by
    changing the database. It has 3 way of working. It can simply mapped reads 
    against an indexed reference database (A), using different mapping tools
    (bowtie2, bwa and novoalign) and different counting methods (best, ex-aequo
    and shared), for non redundant databases. Either perform, in addition, a
    variant calling step, by mapping reads against a clustered redundant database,
    to able an accurate quantification at a gene level (B). It also provide two 
    presets modes (C) for the identification of species (option : --species) 
    and resistance genes (option : --resistance) by mapping reads against the 
    provided databases, RefMG.v1. for bacterial species, ResFinder for resistance 
    genes.
    
    This module performs analysis of metagenomic data (sequencing data of a 
    mixture of bactera), and his adapted for a SGE cluster use. The module 
    operates on both paired-ends and single-ends data, compressed input is 
    supported and auto-detected from the file name (.gz)
    
    Use '%(prog)s --help' to see all command-line options.""",
        
        epilog = color.BOLD + "\n\nEXAMPLES:\n\n" + color.END +
"""   python %(prog)s -h                       to print help

   python %(prog)s mapping -i data              for mapping reads against a
                           -o output            database with bowtie2 
                           -db database_index   count alignments with shared
                           -e email@pasteur.fr 
                           -q hubbioit  
                           --bowtie2 
                           --shared

   python %(prog)s variant -i data              for mapping reads against a 
                           -o output            database and perform variant 
                           -db database_index   calling
                           -fa database.fa 
                           -matrice VCF_matrices/ 
                           -t 4 
                           -e email@pasteur.fr 
                           -q hubbioit 
                           --bowtie2 
                           --shared

   python %(prog)s mode -i data                 To use the preset mode for 
                        -o output               bacterial species quantification
                        -t 8 
                        -e email@pasteur.fr 
                        -q hubbioit 
                        --species

   python %(prog)s mode -i data                 To use the preset mode for
                        -o output               resistance genes quantification
                        -e email@pasteur.fr 
                        -q hubbioit 
                        --resistance 
""")
    
    
    subparsers = parser.add_subparsers(title = 'Working modes')
    
    # Mappping commands
    mapping_parser = subparsers.add_parser('mapping',
        usage= color.BOLD + '\r       ' +
"""
             MBMA - A Mapping Based Microbiome Analysis tool for
                species and resistance genes quantification  

SYNOPSIS - MBMA mapping
""" + color.END + """
    MBMA mapping -i input_dir -o output_dir -q SGE_queue -e email@pasteur.fr 
                 -db database_index [--bowtie2 | --bwa | --novo] 
                 [--best | --ex_aequo | --shared] [options]
    

    options : [--threads NN] [--mode {SE,PE}] [--jobname STR] [--strict] [-h]

""",
        formatter_class = RawTextHelpFormatter,
        help='Map reads against a database and count mapped reads',
        description=color.BOLD + """
DESCRIPTION

    MBMA mapping""" + color.END + """ map sequencing reads against an indexed database.
    It can use different mapping tools (bowtie2, bwa and novoalign) 
    and different counting methods (best, ex-aequo and shared).
    
    This module performs analysis of metagenomic data (sequencing data of a 
    mixture of bactera), and his adapted for a SGE cluster use. The module 
    operates on both paired-ends and single-ends data, compressed input is 
    supported and auto-detected from the file name (.gz)""",
        epilog = color.BOLD + "EXAMPLES:\n\n" + color.END +
"""   MBMA mapping -h                   to print help

   MBMA mapping -i data              for mapping reads against a
                -o output            database with bowtie2 
                -db database_index   count alignments with shared
                -e email@pasteur.fr 
                -q hubbioit  
                --bowtie2 
                --shared""")

    mapping_parser.add_argument("-i", "--indir",
        type = is_dir,
        required = True,
        action = "store",
        metavar = "DIR",
        help = """Path to a directory containing all input FASTQ files.
FASTQ files must have the extension .FASTQ or .FQ.
For paired-ends files must be named : 
    file_R1.[fastq/fq] & file_R2.[fastq/fq] 
                       or
    file_1.[fastq/fq] & file_2.[fastq/fq]\n\n""")
    
    mapping_parser.add_argument("-o", "--outdir",
        type = dir_not_exists,
        required = True,
        metavar = "DIR",
        action = "store",
        help = "Path to the output directory (shouldn't exist before, \n"
               "will be created).\n\n")
    
    mapping_parser.add_argument("-db","--ref_db",
        type = str,
        required = True,
        action = "store",
        metavar = "STRING",
        help = "Prefix of the reference database index for the \n"
               "corresponding mapping tool.\n\n")
    
    mapping_parser.add_argument('-m', "--mode",
        type = is_mode,
        nargs = '?',
        metavar = "PE/SE",
        default = 'PE',
        action="store",
        help = "Reads layout can be single-ends (SE) or paired-ends (PE)\n"
               "  Usage : -m SE or -m PE   [default = 'PE'].\n\n")

    mapping_parser.add_argument("-t", "--threads",
        type = int,
        action = "store",
        nargs = '?',
        default = 1,
        metavar = "INT",
        help = "Number of threads to use.\n"
               "  Usage: -t 5  [default = 1]\n\n")

    mapping_parser.add_argument("-e", "--email",
        type = str,
        required = True,
        metavar="STRING",
        help = "Email for qsub notification. Should be an pasteur email.\n\n")

    mapping_parser.add_argument("-q", "--queue",
        type = str,
        required = True,
        metavar = "STRING",
        help = "Queue name for SGE qsub, there is a specific queue \n"
               "for your team.\n\n")
    
    mapping_parser.add_argument("-j", "--jobname",
        type = str,
        nargs = '?',
        default = "Job",
        metavar = "STRING",
        help = "The name name for this qub job. Jobname will appear \n"
               "in the qstat box.\n"
               "  Usage : -j dataset1 [default 'Job']\n\n")

    mapping_parser.add_argument("-s", "--strict",
        action = 'store_const',
        const = 'strict',
        metavar = "STRING",
        help = "Increase mapping score threashold for a strict mapping.\n"
               "Threshold of Bowtie2 and BWA are increased of 20 X.\n"
               "  Usage : --strict \n\n")
    
    mapping_tool = mapping_parser.add_mutually_exclusive_group(required=True)
    mapping_tool.add_argument('--bowtie2',
        action = 'store_const',
        const = 'True',
        help = "Use Bowtie2 to map read against reference database\n\n")
    mapping_tool.add_argument('--bwa',
        action='store_const',
        const = 'True',
        help = "Use BWA to map read against reference database\n\n")
    mapping_tool.add_argument('--novo', 
        action='store_const',
        const = 'True',
        help = "Use NovoAlign to map read against reference database\n\n")
    
    count_mode = mapping_parser.add_mutually_exclusive_group(required=True)
    count_mode.add_argument('--best', 
        action = 'store_const',
        const = 'True',
        help = "Use the mode 'Best' to count mapped read. 'Best' count\n"
               "only the best alignment (highest alignment score) for\n"
               "each mapped read\n\n")
    count_mode.add_argument('--ex_aequo', 
        action='store_const',
        const = 'True',
        help = "Use the mode 'Ex-aequo' to count mapped read. It count\n"
               "all the best alignments for each mapped read\n\n")
    count_mode.add_argument('--shared', 
        action='store_const',
        const = 'True',
        help = "Use the mode 'Shared' to count mapped read.\n"
               "'Shared' weights the alignments according to the\n"
               "probability that the alignment is the true point of origin\n"
               "of the read")
    
    # Variant commands
    variant_parser = subparsers.add_parser('variant', 
        usage = color.BOLD + '\r       ' +
"""
             MBMA - A Mapping Based Microbiome Analysis tool for
                species and resistance genes quantification  

SYNOPSIS - MBMA variant
""" + color.END + """
    MBMA variant -i input_dir -o output_dir -q SGE_queue -e email@pasteur.fr
                 -db database_index [--bowtie2 | --bwa | --novo] 
                 [--best | --ex_aequo | --shared] [options] -fa database.fasta 
                 -matrice variants_dir [options]


    options : [--threads NN] [--mode {SE,PE}] [--jobname STR] [--strict] [-h]

""",
        formatter_class = RawTextHelpFormatter,
        help='Call variants from mapping and count reads',
        description=color.BOLD + """
DESCRIPTION

    MBMA variant""" + color.END + """ map sequencing reads against an indexed database and performs
    variant calling by using GATK. This mode, allow an accurate quantification 
    of genes for a redundant database. The database should be first clustered 
    using the script make_matrices.py. Then reads are mapped against the provided
    clustered database. The variant profile from the samples are compared to 
    the variant matrices obtained by the clustering, to identify the allele of 
    a gene.    
    
    This module performs analysis of metagenomic data (sequencing data of a 
    mixture of bactera), and his adapted for a SGE cluster use. The module 
    operates on both paired-ends and single-ends data, compressed input is 
    supported and auto-detected from the file name (.gz)""",
    epilog = color.BOLD + "EXAMPLES:\n\n" + color.END +
"""   MBMA variant -h                       to print help

   MBMMA variant -i data              for mapping reads against a 
                 -o output            database and perform variant 
                 -db database_index   calling
                 -fa database.fa 
                 -matrice VCF_matrices/ 
                 -t 4 
                 -e email@pasteur.fr 
                 -q hubbioit 
                 --bowtie2 
                 --shared""")

    variant_parser.add_argument("-i", "--indir",
        type = is_dir,
        required = True,
        action = "store",
        metavar = "DIR",
        help = """Path to a directory containing all input FASTQ files.
FASTQ files must have the extension .FASTQ or .FQ.
For paired-ends files must be named : 
    file_R1.[fastq/fq] & file_R2.[fastq/fq] 
                       or
    file_1.[fastq/fq] & file_2.[fastq/fq]\n\n""")
    
    variant_parser.add_argument("-o", "--outdir",
        type = dir_not_exists,
        required = True,
        metavar = "DIR",
        action = "store",
        help = "Path to the output directory (shouldn't exist before, \n"
               "will be created).\n\n")
    
    variant_parser.add_argument("-db","--ref_db",
        type = str,
        required = True,
        action = "store",
        metavar = "STRING",
        help = "Prefix of the reference database index for the \n"
               "corresponding mapping tool.\n\n")
    
    variant_parser.add_argument("-fa","--ref_fasta",
        type = is_fasta,
        required = True,
        action = "store",
        metavar = "FILE",
        help = "The reference clustered database in FASTA format.\n\n")

    variant_parser.add_argument("-matrice", "--variant_matrice",
        type = is_dir,
        required = True,
        action = "store",
        metavar = "DIR",
        help = "Path to a directory containing all VCF of the variant matrice\n"
               "of the clusters \n\n")

    variant_parser.add_argument('-m', "--mode",
        type = is_mode,
        nargs = '?',
        metavar = "PE/SE",
        default = 'PE',
        action="store",
        help = "Reads layout can be single-ends (SE) or paired-ends (PE)\n"
               "  Usage : -m SE or -m PE   [default = 'PE'].\n\n")

    variant_parser.add_argument("-t", "--threads",
        type = int,
        action = "store",
        nargs = '?',
        default = 1,
        metavar = "INT",
        help = "Number of threads to use.\n"
               "  Usage: -t 5  [default = 1]\n\n")

    variant_parser.add_argument("-e", "--email",
        type = str,
        required = True,
        metavar="STRING",
        help = "Email for qsub notification. Should be an pasteur email.\n\n")

    variant_parser.add_argument("-q", "--queue",
        type = str,
        required = True,
        metavar = "STRING",
        help = "Queue name for SGE qsub, there is a specific queue \n"
               "for your team.\n\n")
    
    variant_parser.add_argument("-j", "--jobname",
        type = str,
        nargs = '?',
        default = "Job",
        metavar = "STRING",
        help = "The name name for this qub job. Jobname will appear \n"
               "in the qstat box.\n"
               "  Usage : -j dataset1 [default 'Job']\n\n")

    variant_parser.add_argument("-s", "--strict",
        action = 'store_const',
        const = 'strict',
        metavar = "STRING",
        help = "Increase mapping score threashold for a strict mapping.\n"
               "Threshold of Bowtie2 and BWA are increased of 20 X.\n"
               "  Usage : --strict \n\n")
    
    mapping_tool = variant_parser.add_mutually_exclusive_group(required=True)
    mapping_tool.add_argument('--bowtie2',
        action = 'store_const',
        const = 'True',
        help = "Use Bowtie2 to map read against reference database\n\n")
    mapping_tool.add_argument('--bwa',
        action='store_const',
        const = 'True',
        help = "Use BWA to map read against reference database\n\n")
    mapping_tool.add_argument('--novo', 
        action='store_const',
        const = 'True',
        help = "Use NovoAlign to map read against reference database\n\n")
    
    count_mode = variant_parser.add_mutually_exclusive_group(required=True)
    count_mode.add_argument('--best', 
        action = 'store_const',
        const = 'True',
        help = "Use the mode 'Best' to count mapped read. 'Best' count\n"
               "only the best alignment (highest alignment score) for\n"
               "each mapped read\n\n")
    count_mode.add_argument('--ex_aequo', 
        action='store_const',
        const = 'True',
        help = "Use the mode 'Ex-aequo' to count mapped read. It count\n"
               "all the best alignments for each mapped read\n\n")
    count_mode.add_argument('--shared', 
        action='store_const',
        const = 'True',
        help = "Use the mode 'Shared' to count mapped read.\n"
               "'Shared' weights the alignments according to the\n"
               "probability that the alignment is the true point of origin\n"
               "of the read")

    # Presets modes commands
    presets_parser = subparsers.add_parser('mode', 
        usage = color.BOLD + '\r       ' +
"""
             MBMA - A Mapping Based Microbiome Analysis tool for
                species and resistance genes quantification  

SYNOPSIS - MBMA mode
""" + color.END + """
    MBMA mode -i input_dir -o output_dir -q SGE_queue -e email@pasteur.fr
              [--species | --resistance]

    options : [--threads NN] [--mode {SE,PE}] [--jobname STR] [--strict] [-h]

""",
        formatter_class = RawTextHelpFormatter,
        help=("Presets modes for bacterial species and resistance \n"
             "genes quantification"),
        description=color.BOLD + """
DESCRIPTION

    MBMA mode""" + color.END + """ identifies and quantifies species (option : --species) 
    and resistance genes (option : --resistance) by mapping reads against the 
    provided databases, RefMG.v1. for bacterial species, ResFinder for resistance 
    genes. Reads are mapped using a higher threshold (--strict) by Bowtie2, 
    then alignments found are counted using Shared. 
    
    This module performs analysis of metagenomic data (sequencing data of a 
    mixture of bactera), and his adapted for a SGE cluster use. The module 
    operates on both paired-ends and single-ends data, compressed input is 
    supported and auto-detected from the file name (.gz)""",
        epilog = color.BOLD + "EXAMPLES:\n\n" + color.END +
"""   MBMA mode -h                                  To print help

   MBMA mode -i data                            To use the preset mode for 
                        -o output               bacterial species quantification
                        -t 8 
                        -e email@pasteur.fr 
                        -q hubbioit 
                        --species

   MBMA mode -i data                            To use the preset mode for
                        -o output               resistance genes quantification
                        -e email@pasteur.fr 
                        -q hubbioit 
                        --resistance 
""")
    
    presets_parser.add_argument("-i", "--indir",
        type = is_dir,
        required = True,
        action = "store",
        metavar = "DIR",
        help = """Path to a directory containing all input FASTQ files.
FASTQ files must have the extension .FASTQ or .FQ.
For paired-ends files must be named : 
    file_R1.[fastq/fq] & file_R2.[fastq/fq] 
                       or
    file_1.[fastq/fq] & file_2.[fastq/fq]\n\n""")
    
    presets_parser.add_argument("-o", "--outdir",
        type = dir_not_exists,
        required = True,
        metavar = "DIR",
        action = "store",
        help = "Path to the output directory (shouldn't exist before, \n"
               "will be created).\n\n")
    
    presets_parser.add_argument('-m', "--mode",
        type = is_mode,
        nargs = '?',
        metavar = "PE/SE",
        default = 'PE',
        action="store",
        help = "Reads layout can be single-ends (SE) or paired-ends (PE)\n"
               "  Usage : -m SE or -m PE   [default = 'PE'].\n\n")

    presets_parser.add_argument("-t", "--threads",
        type = int,
        action = "store",
        nargs = '?',
        default = 1,
        metavar = "INT",
        help = "Number of threads to use.\n"
               "  Usage: -t 5  [default = 1]\n\n")

    presets_parser.add_argument("-e", "--email",
        type = str,
        required = True,
        metavar="STRING",
        help = "Email for qsub notification. Should be an pasteur email.\n\n")

    presets_parser.add_argument("-q", "--queue",
        type = str,
        required = True,
        metavar = "STRING",
        help = "Queue name for SGE qsub, there is a specific queue \n"
               "for your team.\n\n")
    
    presets_parser.add_argument("-j", "--jobname",
        type = str,
        nargs = '?',
        default = "Job",
        metavar = "STRING",
        help = "The name name for this qub job. Jobname will appear \n"
               "in the qstat box.\n"
               "  Usage : -j dataset1 [default 'Job']\n\n")
    
    work_mode = presets_parser.add_mutually_exclusive_group(required=True)
    work_mode.add_argument('--species', 
        action = 'store_const',
        const = 'True',
        help = "Prediction of species using RefMG.v1 database\n")
    
    work_mode.add_argument('--resistance', 
        action='store_const',
        const = 'True',
        help = "Prediction of resistance genes using ResFinder database")
    
    return parser