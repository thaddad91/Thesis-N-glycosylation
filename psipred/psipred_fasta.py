#!/usr/bin/env python3

"""
Author:     Thierry Haddad
Description: Create protein-level FASTA file of the glycosylated and non-glyco-
             sylated proteins, for use with psipred.
"""

import pickle, sys

loc = "/mnt/scratch/hadda003/"

def unpickle_ara_d():
    with open(loc+"negative_set/pos_neg_glycosites", "rb") as f:
        ara_d = pickle.load(f)
    if type(ara_d) == dict:
        create_proteinlevel_fasta(ara_d)

def create_proteinlevel_fasta(ara_d):
    with open("psipred_fasta.fasta", "w") as f:
        for protein in ara_d:
            if len(ara_d[protein]["pos"]) > 0:
                header = ">{}_POS\n".format(protein)
                seq = ara_d[protein]["sequence"]+"\n"
            else:
                header = ">{}_NEG\n".format(protein)
                seq = ara_d[protein]["sequence"]+"\n"
            f.write(header)
            f.write(seq)
    print("Done.")

if __name__ == "__main__":
    unpickle_ara_d()
