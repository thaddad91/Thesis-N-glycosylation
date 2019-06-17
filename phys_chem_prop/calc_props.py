#!/usr/bin/env python3

"""
Author:         Thierry Haddad
Description:    Runs the protein sequences, for neg and pos labelled proteins, through ProtParam's
                Protein Analysis, returning physicochemical properties per protein. These will then
                be used as features for ML.
"""

from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pickle
import sys

loc = "/mnt/scratch/hadda003/"

def unpickle_ara_d():
    """Deserializes the arabidopsis directory with proteins."""
    with open(loc+"negative_set/pos_neg_glycosites", "rb") as f:
        ara_d = pickle.load(f)
    if type(ara_d) == dict:
        return ara_d
    else:
        print("Failed to unpickle ara_d, not a dict!")
        sys.exit(1)

def physchem_props(ara_d):
    """Calculate the physicochemical properties per protein in ara_d."""
    c = 0
    g = 0
    for protein in ara_d:
        seq = ara_d[protein]["sequence"]
        # Calculates the properties
        if "X" in seq:
            continue  # Skip non-usable sequences, only negs
        if '*' in seq:
            if ara_d[protein]["pos"] != []:
                print(protein)
            continue
        a_seq = ProteinAnalysis(seq)
        # Update ara_d with new physchem properties
        results = [
            a_seq.molecular_weight(),
            a_seq.gravy(),
            a_seq.aromaticity(),
            a_seq.instability_index(),
            a_seq.flexibility(),
            a_seq.isoelectric_point(),
            a_seq.secondary_structure_fraction(),
            ]
        keys = [
            "mol_weight",
            "gravy",
            "aromaticity",
            "instab_index",
            "flexi",
            "iso_point",
            "seq_struct",
        ]
        ara_d[protein]["Properties"] = {}
        for k,v in zip(keys, results):
            ara_d[protein]["Properties"][k] = v
    return ara_d

def serialize_ara_d(ara_d):
    """Serialize the arabidopsis dictionary to a file."""
    with open("ara_d_w_properties","wb") as f:
        pickle.dump(ara_d, f)
    return True

def run_wrapper():
    """Wrapper function."""
    ara_d = unpickle_ara_d()
    ara_d = physchem_props(ara_d)
    serialize_ara_d(ara_d)
    k = list(ara_d.keys())[0]
    print(ara_d[k])

if __name__ == "__main__":
    run_wrapper()