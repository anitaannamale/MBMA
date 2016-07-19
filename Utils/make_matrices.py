#! /usr/bin/env python
# -*- coding: utf8 -*-

"""
"""

__author__ = "Anita Annamalé"
__version__  = "2.0"
__copyright__ = "copyleft"
__date__ = "2016/07"

#-------------------------- MODULES IMPORTATION -------------------------------#

import sys
import argparse
import os
import subprocess
import shlex

#-------------------------- FUNCTIONS DEFINITION ------------------------------#

def is_fasta(path):
    """
    Check if a path is an existing FASTA file.

    Args:
        path [STR] = Path to the file

    Returns:
        abs_path [STR] = absolute path
        or quit
    """
    abs_path = os.path.abspath(path) # get absolute path

    # if not a file
    if not os.path.isfile(abs_path):
        # check if it is a path to a dir
        if os.path.isdir(abs_path):
            msg = "{0} is a directory not a file.".format(abs_path)
        # else the path doesn't not exist
        else:
            msg = "The path {0} does not exist.".format(abs_path)
        raise argparse.ArgumentTypeError(msg)

    else :
        ext = os.path.splitext(abs_path)[1] # get the file extension
        if (ext != 'fa') or (ext != 'fasta'):
            msg = ("{0} isn't a fasta file. "
                  "Please, provide a database in FASTA format with extension \
.fasta or .fa".format(abs_path))

    return abs_path
    
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
    with open("{0}/submission.sh".format(param['output']), "wt") as out:
        out.write(temp.format(**param))
    # close the template
    filein.close()
    

#---------------------------- PROGRAMME MAIN ----------------------------------#

if __name__ == '__main__' :
    
    parser = argparse.ArgumentParser(description = "VCF matrices generation \
program ")
    
    parser.add_argument('database', 
                        metavar = 'database',
                        type = is_fasta,
                        help = "input database in fasta format to clusterize\n")
    
    parser.add_argument("-e", "--email",
                        type = str,
                        required = True,
                        metavar="STRING",
                        help = "Email for qsub notification. Should be an \
pasteur email.\n\n")
    
    parser.add_argument("-q", "--queue",
                        type = str,
                        required = True,
                        metavar = "STRING",
                        help = "Queue name for SGE qsub, there is a specific \
queue \nfor your team.\n\n")

    parser.add_argument("-j", "--jobname",
                        type = str,
                        nargs = '?',
                        default = "Job",
                        metavar = "STRING",
                        help = "The name name for this qub job. Jobname will \
appear \nin the qstat box.\n  Usage : -j dataset1 [default 'Job']\n\n")
        
    parser.add_argument('-output',
                        action = 'store',
                        metavar = 'DIR', 
                        type = is_dir,
                        default = os.getcwd(),
                        help = "output path\n")
 
    parser.add_argument('-id_threshold',
                        action = 'store',
                        metavar = 'float', 
                        type = float,
                        default = 0.90,
                        help = "sequence identity threshold for clustering two \
sequence. Must be between 0 and 1 [default : 0.90]\n")
    
    
    
    
    # command-line parsing
    if not len(sys.argv) > 1:
        parser.print_help()
        exit(1)

    try :
        args = parser.parse_args()
        args = vars(args)
        
    except:
        exit(1)
    
    # clean the dictionary
    for key, value in args.items():
            if value == None:
                del args[key]
                
    # get script location
    script_path = os.path.dirname(os.path.abspath(sys.argv[0]))
    
    # complete dictionary
    database = args['database']
    basename_wo_ext = os.path.splitext(os.path.basename(database))[0]
    clusterized_db = '{0}/cdhit_{1}_{2}.fa'.format(args['output'], 
                                                   args['id_threshold'], 
                                                   basename_wo_ext)
    args['clustered_db'] = clusterized_db
    args['script_loc'] = script_path
    args['clstr_fasta'] = '{0}/clstr_fasta'.format(args['output'])
    args['stdout'] = '{0}/{1}.out'.format(args['output'], args['jobname'])
    
    # create a submission script
    submission_template = script_path + '/jobs_make_matrices.sh'
    
    try :
        write_soum(submission_template, args)
    except (IOError, OSError) :
        sys.exit("Template file is missing")
    
    # execute the job
    commandline = "qsub {output}/submission.sh".format(**args)
    arg_split = shlex.split(commandline)
    
    with open("{0}/stdout.txt".format(args['output']),"wt") as out, \
         open("{0}/stderr.txt".format(args['output']),"wt") as err:
        try:
            p = subprocess.call(arg_split, stdout=out, stderr=err)

        except subprocess.CalledProcessError as e:
            p = e.p
            sys.exit(e.returncode)

    errfile = "{output}/stderr.txt".format(**args)
    if os.stat(errfile).st_size != 0:
        with open(errfile, "at") as errfile:
            sys.stderr.write(errfile.read())
        sys.exit(1)    