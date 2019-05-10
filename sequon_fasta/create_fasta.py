#!/usr/bin/env python

"""
Author:     Thierry Haddad
Description: Creates a .FASTA file of all sequons, with their index and class as header
"""

import pickle, sys

loc = "/mnt/scratch/hadda003/negative_set/"

def unpickle():
    """Deserialize dictionary with pos/neg labelled data."""
    with open(loc+"pos_neg_glycosites","rb") as f:
        ara_d = pickle.load(f)
    if type(ara_d) == dict:
        print("Serialization success!")
    else:
        print("Serialization failed!")
        print(type(ara_d))
        sys.exit(1)
    return ara_d

def write_to_fasta(ara_d):
    """Create a .FASTA containing all the sequons"""
    with open("sequon_fasta.fasta", "w") as f:
        for id in ara_d:
            for b in ["pos", "neg"]:
                for seq in ara_d[id][b]:
                    if seq[1] == "" or len(seq[1]) < 15:
                        continue  # Pass empty sequons
                    else:
                        header = ">{}_{}_{}\n".format(id, b, seq[0])
                        f.write(header)
                        f.write(seq[1]+"\n")
    print("Done.")


def run_script():
    """Wrapper function"""
    ara_d = unpickle()
    write_to_fasta(ara_d)

if __name__ == "__main__":
    run_script()