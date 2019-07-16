#!/usr/bin/env python3

"""
Author: Thierry Haddad
Description: Parses the ara_d dictionary. Afterwards, every
             protein will be put through psiblast to generate
             an ASCII PSSM as conservation score feature.
"""

from multiprocessing import Pool
import sys
import subprocess, shlex
import os
import time

loc = '/mnt/scratch/hadda003/individual_fastas/'

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
    # Create command
    #db_loc = "/mnt/nexenta/reference/blast_latest/databases.nobackup/extracted_latest.nobackup/nr"
    db_loc = "/mnt/scratch/hadda003/uniprot_db/uniprot_sprot.fasta"
    #cmd = "psiblast -query {} -db {} -num_threads 5 -num_iterations 2 -out_ascii_pssm {}.txt"  # psiblast
    cmd = "run_psipred.pl -d {} {} -o {}"  # psipred
    cmd = cmd.format(db_loc, loc+fasta, fasta.split('.fa')[0])  # psipred
    #cmd = cmd.format(loc+fasta, db_loc, fasta.split('.fa')[0])  # psiblast
    # Create list to avoid shell=True
    args = shlex.split(cmd)
    # Run the command
    #print("Running ", fasta)
    p = subprocess.check_call(args, shell=False)

def run_wrapper():
    """One-time only for fasta cutting"""
    #fasta = read_fasta()
    #fasta_cutter(fasta)
    # Add all individual fasta files to iterable
    fastas = []
    for file_ in os.listdir(loc):
        if file_.endswith(".fa"):
            fastas.append(str(file_))
    print(len(fastas))
    p = Pool(8)
    p.map(run_psiblast, fastas)

if __name__ == "__main__":
	run_wrapper()
