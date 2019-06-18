#!/usr/bin/env python3

"""
Author:         Thierry Haddad
Description:    Create k-mers for the subset of pos/neg proteins, but does this on protein-level
                instead of sequon-level.
"""

import itertools
import sys
import os

loc = "/mnt/scratch/hadda003/"

def read_table():
    """Reads the subset table and returns it."""
    with open(loc+"sequon_physchem_props.txt","r") as f:
        data = f.readlines()
    return data

def create_kmer_dict():
    """Create dictionary of all possible kmers, set default occurrence to 0"""
    aa = list("ACDEFGHIKLMNPQRSTVWY")  # All 20 natural amino acids
    perm = list(itertools.combinations_with_replacement(aa, 6))  # All triplet permutations
    kmers = []
    for p in perm:
        kmer = "".join(p)
        kmers.append(kmer)
    return kmers

def count_occurrences(data, kmers):
    """For every triplet, count the occurrences in the sequons."""
    lines = []
    header = data[0].rstrip("\n")+"\t"+"\t".join(kmers)+"\n"
    total = len(data[1:])
    current = 0
    for line in data[1:]:
        current += 1
        print(current, "/", total)
        split_line = line.rstrip().split('\t')
        seq = split_line[4]
        # Calculates the properties
        if "X" in seq or '*' in seq:
            continue  # Skip non-usable sequences, only negs
        counts = []
        for kmer in kmers:
            c = str(seq.count(kmer))
            counts.append(c)
        line = line.rstrip("\n")
        line += "\t"+"\t".join(counts)
        lines.append(line)
    with open("subset_kmer_list.txt","w") as kl:
        kl.write(header)
        for line in lines:
            kl.write(line+"\n")

def run_wrapper():
    """Wrapper function."""
    data = read_table()
    kmers = create_kmer_dict()

    count_occurrences(data, kmers)


if __name__ == "__main__":
    run_wrapper()