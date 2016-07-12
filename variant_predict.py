#! /usr/bin/env python
# -*- coding: utf8 -*-

"""
variant_predict.py : Script that 
                        1) Create a matrice merging all the variant matrix from
                           each cluster
                        2) Create a variant profile from a sample VCF file
                        3) Compare variant profile to variant matrice 
                            to identify and quantify genes
                        4) Rewrite a count table
"""

__author__ = "Anita Annamalé"
__version__  = "0.2"
__copyright__ = "copyleft"
__date__ = "2016/07"

#-------------------------- MODULES IMPORTATION -------------------------------#

import sys
from pysam import VariantFile
from collections import defaultdict
import os
import re

#-------------------------- FUNCTIONS DEFINITION ------------------------------#

def parse_vcf(filename):
    """
    Function that parse a database VCF file obtained by a variant calling using
    a multiple alignment file. It parses the VCF file and output a matrix 
    containing all the variant at each snp position of all the clustered genes
    
    Args :
        filename [string] = VCF filename
        
    Returns:
        name [string] = representative gene name 
        index [dict] = a dictionary containing index of snp position in list:
                       key : snp position
                       value : index of the snp in the list of the dict versions
        matrix [dict] = dictionary containing all variations
                          key : clustered gene
                          value : list of the nucleotide variation
    """
    # open VCF file
    vcf = VariantFile(filename)
    # initialise
    index = {}  
    matrix = defaultdict(list)
    i = 0 # index of snp
    name = 0

    for rec in vcf.fetch():
        name = rec.chrom # representative gene name
        # get the snp position (rec.pos) and his index (i)
        index[rec.pos] = i  
        i += 1
        # creation of the matrix of a cluster, gene are the different clustered 
        # genes, obj contain information about the snp
        for gene, obj in rec.samples.items():
            snp = obj.allele_indices[0]
            if snp != -1:
                matrix[gene].append(rec.alleles[snp]) 
            else: # if deletion
                matrix[gene].append('') 
            
    return name, [index, matrix]


def get_variants(filename):
    """
    Function that parse the sample VCF file. This function get snp found in the
    representative genes, and uses the tag 'AD', a list containing the number of
    read mapped for reference and alternative variant.
    
    Args: 
        filename [string] = sample filename
        
    Returns:
        var [dict] = contain snp variation informations of representative genes 
                    key : representative gene name
                    value : variant [dict] containing snp position as key,
                            and a list of (nucleotide variant, aligned reads 
                            number)
    """
    # open VCF file
    vcf = VariantFile(filename)
    # initialise
    var = {}
    flag = 0
    
    for rec in vcf.fetch():
        # only for the first record, set variable name
        if flag == 0:
            name = rec.chrom # rec.chrom is the representative gene name
            variant = defaultdict(list)
            flag = 1
        # if snp are found in another representative gene
        if rec.chrom != name:
            var[name] = variant # store the variant
            name = rec.chrom # change the representative gene name
            variant = defaultdict(list) # create a new variant dictionnary
        # read the snp informations
        for gene, obj in rec.samples.items():
            i = 0
            if 'AD' in obj:
                for nb in obj['AD']:
                    if nb != 0:
                        variant[rec.pos].append((rec.alleles[i], nb))
                    i +=1
    
    return var



def predict(database, sample_var):
    """
    Function that identifies and quantifies genes.

    Args: 
        database [dict] = contain all genes and their variations
        sample_var [dict] = contain sample snp variations

    Returns:
        res [dict] = Containas key representative, genes and their percentage as
        			 value.
    """
    ## STEP 1 : GET only representative genes whose variant and are present
    ##              in the database --------------------------------------------
    sample_represent = sample_var.keys() # representative presenting variation
    db_represent = database.keys() 
    var_representatives = set(db_represent).intersection(sample_represent)
    
    # initialise
    res = {}
    qt = {}
    
    for represent in var_representatives:
        ref = database[represent] # get the representative variant matrix
        ref_index = ref[0]
        ref_matrix = ref[1]
        positions = ref_index.keys() # all snp positions 
        var_pos = sample_var[represent].keys() # all snp position in sample
        commun_pos = set(positions).intersection(var_pos) # commun snp position
        versions = ref_matrix.keys() # all clustered genes name
        
        ## STEP 2 : Identification of versions of a representative -------------
        for pos in positions:
            idx = ref_index[pos] # get the snp index
            if pos not in var_pos:
                snp = ref_matrix[represent][idx] # variant
                versions = [gene for gene in versions 
                            if ref_matrix[gene][idx] == snp]
            else:
                snp = [ snp[0] for snp in sample_var[represent][pos] ] # variant
                versions = [gene for gene in versions 
                            if any(ref_matrix[gene][idx] == variant 
                                for variant in snp)]
                if len(versions) == 0:
                    versions = [represent]

        ## STEP 3 : Quantification of versions --------------------------------- 
        if len(versions) == 1:
            if represent != versions[0]:
                res[represent] = versions[0]
        else:
            for c_pos in commun_pos:
                idx = ref_index[c_pos]
                snp = [ snp for snp in sample_var[represent][c_pos] ]
                # get total number of mapped reads at that position
                total = sum([val[1] for val in snp])
                # quantification
                for variant, count in snp:
                    quant_versions = [ gene for gene in versions 
                                      if ref_matrix[gene][idx] == variant ]
                    if len(quant_versions) == 1:
                        qt[quant_versions[0]] = round(count/float(total),3)
            res[represent] = qt
    
    return res


def write_count_table(filename, prediction, fileout):
    """
    Function that take the old count table, and the quantification result to 
    write a new count table.
    
    Args: 
        filename [str] = current count table filename
        prediction [dict] = contain new prediction results
        fileout [str] = name of the new count table.

    No returns
    """
    with open(fileout, "wt") as fileout:
        with open(filename, "rt") as filein:
            for line in filein:
                columns = line.split()
                gene = columns[0]
                if gene in prediction.keys(): 
                    new_gene_info = prediction[gene]
                    # if representative gene have been replaced by one gene
                    if type(new_gene_info) == str:
                        fileout.write('{0}\t{1}\t{2}\n'.format(new_gene_info,
                                                               columns[1], 
                                                               columns[2]))
                    # if representative gene have been replaced by two or more 
                    # gene
                    elif type(new_gene_info) == dict:
                        total = float(columns[2])
                        for name, perc in new_gene_info.items():
                            new_count = int(total*perc)
                            fileout.write('{0}\t{1}\t{2}\n'.format(name,
                                                                   columns[1],
                                                                   new_count))
                else :
                    fileout.write(line)
            

#--------------------------------- MAIN ---------------------------------------#

if __name__ == '__main__' :
    # command-line parsing
    if not len(sys.argv) == 5:
        print ("usage: python database.py <dossier_vcf> <fichier.vcf> "
               "<count_table.txt> <new_table>.txt")
        exit(1)
    dossier = os.path.abspath(sys.argv[1])
    sample_vcf = os.path.abspath(sys.argv[2])
    table = os.path.abspath(sys.argv[3])
    out_table = os.path.abspath(sys.argv[4])

    ## STEP 1 : Database creation ----------------------------------------------
    db = {}
    # get all files in the variant matrices directory 
    onlyfiles = [os.path.join(dossier, f) for f in os.listdir(dossier) 
                 if os.path.isfile(os.path.join(dossier, f))]
    # Parse all the VCF from the directory to construct the database
    for vcf in onlyfiles:
        gene, info = parse_vcf(vcf)
        db[gene] = info

    ## STEP 2 ; Parse sample VCF file ------------------------------------------
    variation = get_variants(sample_vcf)

    ## STEP 3 : Identify and Quantify genes ------------------------------------
    results  = predict(db, variation)

    ### STEP 4 : Write the new count table -------------------------------------
    write_count_table(table, results, out_table)