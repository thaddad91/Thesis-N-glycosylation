#!/usr/bin/env python3

"""
Author: Thierry Haddad
Description: Parses the ara_d dictionary. Afterwards, every
             protein will be put through psiblast to generate
             an ASCII PSSM as conservation score feature.
"""

from multiprocessing.dummy import Pool as ThreadPool
import sys
import subprocess, shlex

loc = '/mnt/scratch/hadda003/'

def read_fasta():
    """Reads complete FASTA set of all proteins"""
    with open(loc+"complete_fasta.fa", "r") as f:
        fasta = f.readlines()
    return fasta

def fasta_cutter(fasta):
    """Cuts complete FASTA into FASTA's with a single entry"""
    for row, line in enumerate(fasta, 1):
        if line.startswith(">"):
            header = line
            file_name = header.strip().lstrip(">")
            sequence = fasta[row]
            with open(loc+"individual_fastas/"+file_name+".fa", "w") as f:
                f.write(header)
                f.write(sequence)

def run_psiblast(fasta, db_loc):
    """Run threaded psiblast for every individual protein FASTA"""
    # Create command
    # TO-DO fasta entry, output
    cmd = "run_psipred.pl -d {} {}".format(db_loc, fasta)
    # Create list to avoid shell=True
    args = shlex.split(cmd)
    # Run the command
    p = subprocess.Popen(args, stdout=subprocess.PIPE, shell=False)

    return answer

def run_wrapper():
    fasta = read_fasta()
    fasta_cutter(fasta)
    sys.exit(1)
    db_loc = "/mnt/nexenta/reference/blast_latest/databases.nobackup/extracted_latest.nobackup/nr"
    pool = ThreadPool(4)
    results = pool.map(run_psiblast(fasta, db_loc), fasta)

if __name__ == "__main__":
	run_wrapper()
