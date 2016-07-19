#! /usr/bin/env python
# -*- coding: utf8 -*-


__author__ = "Anita Annamalé"
__version__  = "1.0"
__copyright__ = "copyleft"
__date__ = "2016/05"

#-------------------------- MODULES IMPORTATION -------------------------------#

import sys
import os

#-------------------------- FUNCTIONS DEFINITION ------------------------------#

def parse_cltr(filename):
    """
    Function that parse a cluster information file from CDHIT-est
    
    Args : 
        filename [STR] = cluster informations file name
    
    Returns:
        annot [DICT] = contains gene as key and his cluster name as value
        clstr_vs_ref [DICT] = contains cluster name as key and cluster 
                              representative gene as value
    """
    # initialize
    annot = {}
    clstr_vs_ref = {}
    
    with open(filename, "rt") as infile:
        for line in infile:
            if line.startswith('>Cluster '):
                clst_nb = line.rstrip()
                clst_nb = clst_nb.replace(" ", "") # number of the cluster
            else :
                genes_info = line.split(" ")
                # annotate each gene to a cluster
                annot[genes_info[1][:-3]] = clst_nb[1:]
                # if gene is the representative gene of the cluster
                if genes_info[2] == '*\n':
                    # annotate each cluster to his reference sequence
                    clstr_vs_ref[clst_nb[1:]] = genes_info[1][1:-3]
                    
    return annot, clstr_vs_ref


def create_multifasta(fasta_file, annot, clstr_vs_ref, outdir):
    """
    Function that create multifasta file in th provided directory (outdir).
    Multifasta files are named by the clusters representatives name, and 
    contains all the genes that are similar to the representatives.
    
    Args:
        fasta_file [STR] = fasta file of the un clusterized database
        annot [DICT] = contains gene as key and his cluster name as value
        clstr_vs_ref [DICT] = contains cluster name as key and cluster 
                              representative gene as value
        outdir [STR] = output directory name
    
    No return
    """
    # initialize
    fasta_data = {}
    
    # open the database fasta file and get all information
    with open(fasta_file, "rt") as multifasta:
        for line in multifasta:
            if line.startswith('>'):
                name = line.rstrip()  # get gene name
                clstr_name = annot[name] # get gene cluster name
                ref_name = clstr_vs_ref[clstr_name]  # get cluster 
                                                     # representative gene name
            # append to the fasta dictionnary 
            fasta_data.setdefault(ref_name, []).append(line)
                                                            
    # write all the multifasta files
    for element in fasta_data: # for representative
        with open('{0}/{1}.fa'.format(outdir,element), "wt") as output:
            output.write("".join(fasta_data[element])) # write sequences
    

#--------------------------------- MAIN ---------------------------------------#

if __name__ == '__main__' :
    
    # command-line parsing
    if not len(sys.argv) == 4:
        msg = "usage: python cluster2multifasta.py <file.clstr> <file.fasta> \
<fasta_dir>"
        print msg 
        exit(1)
        
    # get absolute path
    clstr_name = os.path.abspath(sys.argv[1])
    fasta_name = os.path.abspath(sys.argv[2])
    dir_name = os.path.abspath(sys.argv[3])
    
    # parse the cluster information file from CDHIT
    annotate_gene, get_ref = parse_cltr(clstr_name)
    
    # create a repertory in the working directory
    os.mkdir(dir_name)
    
    # create multifasta files
    create_multifasta(fasta_name, annotate_gene, get_ref, dir_name)
    
    
    
    
