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
    """Create protein-level subset for positive and negatives."""
    lines = []
    header = "ID\tClass\tSequence\n"
    lines.append(header)
    for id in ara_d:
        if len(ara_d[id]["pos"]) > 0:
            line = "{}\t{}\t{}\n".format(id, "pos", ara_d[id]["sequence"])
        else:
            line = "{}\t{}\t{}\n".format(id, "neg", ara_d[id]["sequence"])
        lines.append(line)
    return lines

def write_table(lines):
    """Writes protein-level set to tab-delimited file."""
    with open(loc+"protein_level_classes.txt","w") as f:
        for line in lines:
            f.write(line)

def run_wrapper():
    ara_d = read_pickle()
    lines = create_subset(ara_d)
    write_table(lines)

if __name__ == "__main__":
    run_wrapper()