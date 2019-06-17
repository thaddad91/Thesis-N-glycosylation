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
    with open(loc+"subset_table.txt","r") as f:
        data = f.readlines()
    return data

def physchem_props(data):
    """Calculate the physicochemical properties per protein in ara_d."""
    new_table = []
    header = "ID\tclass\tindex\tsequon\tsequence\tmol_weight\tgravy\taromaticity\tinstab_index\tiso_point\n"
    new_table.append(header)
    for line in data:
        split_line = line.rstrip().split('\t')
        seq = split_line[-2]  # Sequon, not sequence
        # Calculates the properties
        if "X" in seq or '*' in seq or seq == '':
            continue  # Skip non-usable sequences, only negs
        try:
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
        except:
            print(split_line)
            sys.exit(1)
        new_line = line.rstrip()+"\t{}\t{}\t{}\t{}\t{}\n".format(*results)
        new_table.append(new_line)
    return new_table

def write_table(table):
    """Writes nested list to tab-delimited table file."""
    with open(loc+"sequon_physchem_props.txt","w") as f:
        for line in table:
            f.write(line)

def run_wrapper():
    """Wrapper function."""
    data = read_table()
    table = physchem_props(data)
    write_table(table)


if __name__ == "__main__":
    run_wrapper()
