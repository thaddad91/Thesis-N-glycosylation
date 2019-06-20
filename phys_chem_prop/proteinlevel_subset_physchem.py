#!/usr/bin/env python3

"""
Author:         Thierry Haddad
Description:    Runs the protein sequences, for neg and pos labelled proteins, through ProtParam's
                Protein Analysis, returning physicochemical properties per protein. These will then
                be used as features for ML.
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import sys

loc = "/mnt/scratch/hadda003/"

def read_table():
    """Reads the subset table and returns it."""
    with open(loc+"subset_protein_level_table.txt","r") as f:
        data = f.readlines()
    return data

def physchem_props(data):
    """Calculate the physicochemical properties per protein in subset."""
    new_table = []
    header = "id\tclass\tsequence\tmol_weight\tgravy\taromaticity\tinstab_index\tiso_point\n"
    new_table.append(header)
    for line in data[1:]:
        split_line = line.rstrip().split('\t')
        seq = split_line[-1]
        # Calculates the properties
        if "X" in seq or '*' in seq or 'e' in seq:
            print(split_line)
            sys.exit(1)
            continue  # Skip non-usable sequences, only negs

        a_seq = ProteinAnalysis(seq)
        # Update ara_d with new physchem properties
        results = [
            a_seq.molecular_weight(),
            a_seq.gravy(),
            a_seq.aromaticity(),
            a_seq.instability_index(),
            #a_seq.flexibility(),
            a_seq.isoelectric_point(),
            #a_seq.secondary_structure_fraction(),
            ]
        new_line = line.rstrip()+"\t{}\t{}\t{}\t{}\t{}\n".format(*results)
        new_table.append(new_line)
    return new_table

def write_table(table):
    """Writes nested list to tab-delimited table file."""
    with open(loc+"protein_level_physchem_props.txt","w") as f:
        for line in table:
            f.write(line)

def run_wrapper():
    """Wrapper function."""
    data = read_table()
    table = physchem_props(data)
    write_table(table)


if __name__ == "__main__":
    run_wrapper()