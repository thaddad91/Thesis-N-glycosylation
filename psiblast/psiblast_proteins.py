#!/usr/bin/env python3

"""
Author: Thierry Haddad
Description: Parses the ara_d dictionary. Afterwards, every
             protein will be put through psiblast to generate
             an ASCII PSSM as conservation score feature.
"""

from multiprocessing.dummy import Pool as ThreadPool
import sys

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

def run_psiblast(fasta):
    """Run threaded psiblast for every individual protein FASTA"""
    pass

def run_wrapper():
    fasta = read_fasta()
    fasta_cutter(fasta)
    sys.exit(1)
    pool = ThreadPool(4)
    results = pool.map(run_psiblast, fasta)

if __name__ == "__main__":
	run_wrapper()
