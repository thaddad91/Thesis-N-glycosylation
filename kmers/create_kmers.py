#!/usr/bin/env python

"""
Author:      Thierry Haddad
Description: Script to extract k-mers occurrence per sequon for
             unsupervised learning purposes.
"""

import itertools
import pickle
import sys
import os

def create_kmer_dict():
    """Create dictionary of all possible kmers, set default occurrence to 0"""
    aa = list("ACDEFGHIKLMNPQRSTVWY")  # All 20 natural amino acids
    perm = list(itertools.combinations_with_replacement(aa, 3))  # All triplet permutations
    kmers = []
    for p in perm:
        kmer = "".join(p)
        kmers.append(kmer)
    return kmers

def read_pickle():
    """Open and read serialized dictionary containing glycosites, indixes and sequences."""
    with open("pos_neg_glycosites","rb") as f:
        ara_d = pickle.load(f)
    if type(ara_d) == dict:
        print("Serialization succes!")
    else:
        print("Serialization failed!")
        print(type(ara_d))
        sys.exit(1)
    return ara_d

def count_occurrences(ara_d, kmers):
    """For every triplet, count the occurrences in the sequons."""
    lines = []
    header = "glycosite\t"+"\t".join(kmers)+"\n"
    for id in ara_d:
        for b in ["pos", "neg"]:
            for seq in ara_d[id][b]:
                seq = seq[1]  # Sequon
                line = b+"\t"
                counts = []
                for kmer in kmers:
                    c = str(seq.count(kmer))
                    counts.append(c)
                line += "\t".join(counts)
                lines.append(line)
    with open("kmer_list.txt","w") as kl:
        kl.write(header)
        for line in lines:
            kl.write(line+"\n")

def seq_list(ara_d):
    """Create a list of sequons and a positive/negative label.
       Currently unused.
    """
    results = []
    for id in ara_d:
        for b in ["pos", "neg"]:
            for seq in ara_d[id][b]:
                if seq[1] != '':
                    results.append([seq[1], b])  # Append sequon + pos/neg label
    with open("seq_list.txt","w") as sl:
        for result in results:
            sl.write(result[0]+"\t"+result[1]+"\n")

def run_script():
    kmers = create_kmer_dict()
    ara_d = read_pickle()
    #seq_list(ara_d)
    count_occurrences(ara_d, kmers)

if __name__ == "__main__":
    run_script()
