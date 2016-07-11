#! /usr/bin/env python
# -*- coding: utf8 -*-


__author__ = "Anita Annamalé"
__version__  = "1.0"
__copyright__ = "copyleft"
__date__ = "2016/02"

#-------------------------- MODULES IMPORTATION -------------------------------#

import sys
import re
import os
import pysam
import argparse
from collections import defaultdict
import itertools


#-------------------------- FUNCTIONS DEFINITION ------------------------------#

def is_sam(path):
    """
    Check if a path is an existing file.

    Args:
        path [string] = Path to the file

    Returns:
        abs_path [string] = absolute path
        or quit
    """
    abs_path = os.path.abspath(path) # get absolute path

    # if not a directory
    if not os.path.isfile(abs_path):
        # check if it is a path to a file
        if os.path.isdir(abs_path):
            msg = "{0} is a directory not a file.".format(abs_path)
        # else the path doesn't not exist
        else:
            msg = "The path {0} does not exist.".format(abs_path)
        raise argparse.ArgumentTypeError(msg)

    else :
        ext = os.path.splitext(abs_path)[1]
        if ext != 'sam':
            msg = ("{0} isn't a sam file. "
                  "Please, provide a file with extension .sam".format(abs_path))

    return abs_path


def path_not_exists(path):
    """
    Check if a path is not exist.

    Args:
        path [string] = Path to check

    Returns:
        abs_path [string] = absolute path if it doesn't exist
        or quit
    """
    abs_path = os.path.abspath(path)

    # if path exists
    if os.path.exists(abs_path):
        msg = "The path {0} already exist. "
        "Please give another file name".format(abs_path)
        raise argparse.ArgumentTypeError(msg)

    return abs_path


def create_bam(filename):
    bamfile = os.path.dirname(filename)[:-3] + "bam/" + os.path.basename(filename)[:-3] + "bam" 
    # convert sam to bam
    pysam.view('-Sb',filename, '-o', bamfile, catch_stdout=False)

    return bamfile


def sort_bam(filename):
    sortfile = '{0}/sorted_{1}'.format(os.path.dirname(filename),
                                       os.path.basename(filename)[:-4])
    # sort the bam file
    pysam.sort(filename, sortfile)
    
    return sortfile + ".bam"

    
def only_mapped(filename):
    mappedfile = '{0}/filtered_{1}'.format(os.path.dirname(filename),
                                       os.path.basename(filename)[7:])
    # get only mapped reads
    pysam.view('-b', '-F', '4', filename, '-o', mappedfile, catch_stdout=False)
    
    return mappedfile


def write_table(filename, outfile):
    
    # index the bam file
    pysam.index(filename)
    # create count table
    table = pysam.idxstats(filename)
    
    with open(outfile, 'wt') as out:
        for line in table:
            out.write(line)


def read_bam(filename, count, bam=False):
    tmp_score = dict()
    tmp_genomes = defaultdict(list)
    samfile = pysam.AlignmentFile(filename, "rb")
    
    if (count == "ex-aequo") or bam:
        reads = defaultdict(list)
        new_filename = os.path.dirname(filename) + "/unsorted_filtered_" + os.path.basename(filename)
    if count == "shared":
        references = samfile.references
        lengths = samfile.lengths
        database = dict(zip(references, lengths))
    
    for element in samfile:
        if not element.is_unmapped:
            if element.has_tag('AS'):
                if element.is_read1:
                    read_id = '{0}_1'.format(element.qname)
                else :
                    read_id = '{0}_2'.format(element.qname)
                score = element.get_tag('AS')
                prev_score = tmp_score.get(read_id, score)
                if prev_score == score:
                    tmp_score[read_id] = score
                    if element.reference_name not in tmp_genomes[read_id]:
                        tmp_genomes[read_id].append(element.reference_name)
                        if (count == "ex-aequo") or bam:
                            reads[read_id].append(element)
                elif prev_score < score:
                    tmp_score[read_id] = score
                    tmp_genomes[read_id] = [element.reference_name]
                    if (count == "ex-aequo") or bam:
                        reads[read_id] = [element]
            else :
                sys.exit("[FATAL error] Parsing sam file : no optional 'AS' field.\n")

    samfile.close()

    if (count == "ex-aequo") or bam:
        read_list = list(itertools.chain(reads.values()))
        merged_list = list(itertools.chain.from_iterable(read_list))

        ex_aequo_reads = pysam.AlignmentFile(new_filename, "wb", template=samfile)
        for element in merged_list:
            ex_aequo_reads.write(element)
        ex_aequo_reads.close()
        if count == "ex-aequo":
            return new_filename
        if bam:
            sortfile = os.path.dirname(filename) + "/filtered_" + os.path.basename(filename)[:-4]
            pysam.sort(new_filename, sortfile)
        
    if (count =="shared") or bam:
        for (key,values) in tmp_genomes.iteritems():
            tmp_genomes[key] = list(set(values))   
        return database, tmp_genomes


def uniq_from_mult(genome_dict, unique_dict):
    """
    Function that filter unique reads from all reads. Multiple reads are 
    reads that map to more than one genome. And Unique reads are reads that map 
    only on one genome.

    Args:
        genome_dict [dict] = dictionary containing reads as key and a list of 
                             genome where read mapped as value
        unique_dict [dict] = contains for each reference genome the number of 
                             unique reads

    Returns:
        genome_dict [dict] = the dictionary without unique reads
        unique_dict [dict] = nb of unique read of each reference genome
    """
    unique_reads = []

    for key in genome_dict:
        # if read mapp with best score to only one genome:
        if (len(genome_dict[key]) == 1):
            genome = ''.join(genome_dict[key]) # get genome id
            unique_dict[genome] += 1 # add to unique read dictionnary
            unique_reads.append(key)
    
    for key in unique_reads:
        del genome_dict[key] # delete the key from TMP dictionnary

    return genome_dict, unique_dict


def calculate_Co(genome_dict, unique_dict):
    """
    Calculate genome specific coefficient "Co" for each multiple read.

    Args:
        genome_dict [dict] = Contains for each multiple read, all the reference
                             genome name where he mapped
        unique_dict [dict] = nb of unique read of each reference genome

    Returns:
        Co [dict] = contains coefficient values for each couple (read, genome)
    """
    Co = {}
    read_dict = dict()

    for key in genome_dict:
        s = [ unique_dict[genome] for genome in genome_dict[key] ]
        som = reduce(lambda x,y : x+y, s)
        if (som != 0):
            for genome in genome_dict[key]:
                nb_unique = unique_dict[genome]
                if(nb_unique != 0):
                    Co[(key, genome)] = nb_unique/float(som)
                    read_dict.setdefault(genome, []).append(key)

    for genome, reads in read_dict.iteritems():
        read_dict[genome] = list(set(reads))

    return read_dict, Co


def calcul_AbM(unique_dict, read_dict, Co_dict, multiple_dict):
    """
    Calculates multiple reads abundance for each genome

    Args:
        unique_dict [dict] = nb of unique read of each reference genome
        Co_dict [dict] = contains coefficient values for each couple 
                         (read, genome)
        multiple_dict [dict] = abundance of multiple reads for each genome

    Returns:
        multiple_dict [dict] = abundance of multiple reads for each genome
    """
    for genome in read_dict:
        s = [Co_dict[(read, genome)] for read in read_dict[genome]]
        som = sum(s)
        multiple_dict[genome] = som

    return multiple_dict


def calcul_AbS(unique_dict, multiple_dict):
    """
    Calculates the abundance of a each genome

    Args:
        unique_dict [dict] = nb of unique read of each reference genome
        multiple_dict [dict] = abundance of multiple reads for each genome

    Returns:
        abundance_dict [dict] = contains abundance of each reference genome 
    """
    abundance_dict = {}

    for genome in unique_dict:
        abundance_dict[genome] = unique_dict[genome] + multiple_dict[genome]

    return abundance_dict


def write_stat(output, abundance_dict, database):
    """
    Write count table
    
    Args:
        output [string] = output filename
        abundance_dict [dict] = contains abundance of each reference genome
        database [dict] = contrains length of each reference genome
    """
    with open(output, 'wt') as out:
        for genome, abundance in abundance_dict.items():
            out.write('{0}\t{1}\t{2}\n'.format(genome, database[genome], abundance))



#---------------------------- PROGRAMME MAIN ----------------------------------#

if __name__ == '__main__' :

    parser = argparse.ArgumentParser(description = " Reads counting Program ")
    
    parser.add_argument('samfile', 
                        metavar = 'SAMFILE',
                        type = is_sam,
                        help = "input bam file to parse\n")
    
    parser.add_argument('output', 
                        metavar = 'OUTPUT',
                        type = path_not_exists,
                        help = "Path to output the count table\n")
    
    parser.add_argument('-bam',
                        action = 'store_const',
                        const = 'bam',
                        help = "Write the BAM file\n")

    counting_method = parser.add_mutually_exclusive_group(required=True)
    
    counting_method.add_argument('--best',
                                 action = 'store_const',
                                 const = 'best',
                                 help = "Count each reads only once. "
                                         "If a read map to multiple location with"
                                         "a best score, one location is chosen "
                                         "randomly\n")
    
    counting_method.add_argument('--ex_aequo',
                                 action='store_const',
                                 const = 'ex_aequo',
                                 help = "Count all best hits. "
                                          "If a read have multiple best hits, "
                                          "all best hits are taken into "
                                          "account\n")

    counting_method.add_argument('--shared',
                                 action='store_const',
                                 const = 'shared',
                                 help = "Weights each reads according to the"
                                         "probability that the alignment does "
                                         "not correspond to the read's true "
                                         "point of origin\n")
    
    # command-line parsing
    if not len(sys.argv) > 2:
        parser.print_help()
        exit(1)

    try :
        args = parser.parse_args()
        args = vars(args)
        
    except:
        exit(1)
    
    
    # PARSING PARAMETERS -------------------------------------------------------
    
        # input & output
    samfile = args['samfile']
    output = args['output']
        
        # count mode
    
    # clean the dictionary
    for key, value in args.items():
            if value==None:
                del args[key]

    # First Step: Create BAM file
    bamfile = create_bam(samfile)
    
    # Second Step : Count reads & create count table
    if 'best' in args:
        sort_file = sort_bam(bamfile)            # sort bamfile
        mapped_file = only_mapped(sort_file)     # mapped reads only
        write_table(mapped_file, output)        # write count table
    
    elif 'ex_aequo' in args:
        filtered_bam = read_bam(bamfile, "ex-aequo") # create a bamfile
        sort_file = sort_bam(filtered_bam)
        write_table(sort_file, output)

    else : # shared
        # parsing du fichier bam
        if 'bam' in args:
            bam=True
        else:
            bam = False
        
        db, genomes = read_bam(bamfile, "shared", bam)     

        # create dictionaries
        unique = dict.fromkeys(db, 0)
        multiple = dict.fromkeys(db, 0)
        
        # Filter unique reads from multiple reads
        genomes, unique = uniq_from_mult(genomes, unique)
        
        # For multiple reads
        reads, coef_read = calculate_Co(genomes, unique)   # calculate Co 
        multiple = calcul_AbM(unique, reads, coef_read, multiple) # calculate 
                                                                   # abundance
         
        # Calculate reference abundance & write count table
        abundance = calcul_AbS(unique, multiple)
        write_stat(output, abundance, db)
