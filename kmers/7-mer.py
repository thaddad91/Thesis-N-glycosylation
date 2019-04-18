#!/usr/bin/env python3

"""
Author:      Thierry Haddad
Description: Script to retrieve all k-mers (k=7) from the glycosite dataset. 
             However, unlike the 3-mers these will not just shift in index 
             location, but rather represent the 7 residues next to the
             Asparagine in the Nx[S/T] motif.
"""

import pickle

# Scratch folder where data files are located
data_dir = "/mnt/scratch/hadda003/kmers/"

def read_pickle():
    """Read serialized bytestream of glycosite table."""
    with open(data_dir+"pos_neg_glycosites","rb") as f:
        ara_d = pickle.load(f)
    if type(ara_d) == dict:
        print("Serialization succes!")
    else:
        print("Serialization failed!")
        print(type(ara_d))
        sys.exit(1)
    return ara_d

def count_occurrences(ara_d):
    """Count all occurrences of 7-mers in the 15 residue window."""
    kmers = {}  # Dict of all found 7-mers
    for id in ara_d:
        bool = ["pos", "neg"]
        for b in bool:
            for seq in ara_d[id][b]:
                for kmer in [seq[1][:7], seq[1][8:]]:
                    if len(kmer) < 7:  # Pass too short k-mers (edges)
                        continue
                    if kmer not in kmers:
                        kmers[kmer] = {"pos":0, "neg":0}
                    kmers[kmer][b] += 1
    header = "kmer\tpos\tneg\n"

    with open("7mer_list.txt","w") as kl:
        kl.write(header)
        for kmer in kmers:
            kl.write(kmer+"\t"+str(kmers[kmer]["pos"])+"\t"+str(kmers[kmer]["neg"])+"\n")

def run_script():
    ara_d = read_pickle()
    count_occurrences(ara_d)

if __name__ == "__main__":
    run_script()
