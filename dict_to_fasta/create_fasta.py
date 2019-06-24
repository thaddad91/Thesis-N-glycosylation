#!/usr/bin/env python

"""
Author:      Thierry Haddad
Description: Convert ara_d dict to .FASTA format.
"""

import itertools
import pickle
import sys
import os

loc = "/mnt/scratch/hadda003/"

def read_pickle():
    """Open and read serialized dictionary containing glycosites, indixes and sequences."""
    with open(loc+"negative_set/pos_neg_glycosites","rb") as f:
        ara_d = pickle.load(f)
    if type(ara_d) == dict:
        print("Serialization succes!")
    else:
        print("Serialization failed!")
        print(type(ara_d))
        sys.exit(1)
    return ara_d

def create_subset(ara_d):
    """Create FASTA file with class label in header."""
    lines = []
    for id_ in ara_d:
        if len(ara_d[id_]["pos"]) > 0:
            line = ">{}_{}\n{}\n".format("pos", id_, ara_d[id_]["sequence"])
        else:
            line = ">{}_{}\n{}\n".format("neg", id_, ara_d[id_]["sequence"])
        lines.append(line)
    return lines

def write_fasta(lines):
    """Write list with headers/sequences to .FASTA file format."""
    with open(loc+"complete_fasta.fa","w") as f:
        for line in lines:
            f.write(line)

def run_wrapper():
    ara_d = read_pickle()
    lines = create_subset(ara_d)
    write_fasta(lines)

if __name__ == "__main__":
    run_wrapper()