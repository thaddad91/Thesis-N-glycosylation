#!/usr/bin/env python3

"""
Author:      Thierry Haddad
Description: Creates a .FASTA file for the ±4400 subset table, destined for
             psi-blast PSSM matrices.
"""

loc = "/mnt/scratch/hadda003/"

def read_table():
    """Write tab-del file of ±4400 subset."""
    with open(loc+"subset_table.txt") as f:
        data = f.readlines()
    write_fasta(data)

def write_fasta(data):
    """Write listed data to .FASTA file."""
    pos = 0
    neg = 0
    with open("subset_fasta.fasta", "w") as f:
        for line in data:
            line = line.split("\t")
            seq = line[3].strip()
            if len(seq) == 15:
                header = ">" + "_".join(line[0:3]) + "\n"
                if "pos" in header:
                    pos += 1
                else:
                    neg += 1
                f.write(header)
                f.write(seq+"\n")
    print("Added {} pos and {} neg.".format(pos, neg))

if __name__ == "__main__":
    read_table()
