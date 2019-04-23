#!/usr/bin/env python

"""
Author:     Thierry Haddad
Description: Script to calculate the distance between two A.A. sequences using
             a substitution matrix. Afterwards, assign the class of the lowest
             distance to the new sequence, and calculate accuracy.
"""

from Bio.SubsMat.MatrixInfo import blosum62
from Bio import pairwise2
import pickle
import sys

loc = "/mnt/scratch/hadda003/negative_set/"

def unpickle():
    """Deserialize dictionary with pos/neg labelled data."""
    with open(loc+"pos_neg_glycosites","rb") as f:
        ara_d = pickle.load(f)
    if type(ara_d) == dict:
        print("Serialization succes!")
    else:
        print("Serialization failed!")
        print(type(ara_d))
        sys.exit(1)
    return ara_d

def partition_data(ara_d):
    """Create data partitions for train/test sets."""
    pos_set = []
    neg_set = []
    for id in ara_d:
        for b in ["pos", "neg"]:
            for seq in ara_d[id][b]:
                # Append sequon to positive or negative set
                if b == "pos":
                    pos_set.append(seq[1])
                else:
                    neg_set.append(seq[1])
    pos_set = pos_set[:100]  # Debugging
    sub_neg = neg_set[:100]#2200]  # Subset for debugging purposes

    # Train (75%) and test (25%) sets
    train = {
        "pos": pos_set[:int(len(pos_set) * 0.75)],
        "neg": sub_neg[:int(len(sub_neg) * 0.75)],
    }
    test = {
        "pos": pos_set[int(len(pos_set) * 0.75):],
        "neg": sub_neg[int(len(sub_neg) * 0.75):],
    }
    return train, test

def calc_distance(train, test):
    """Calculate distance between seq from test and sequences from train.
    Assign class of lowest distance train sequence to test sequence.
    """
    good = 0
    bad = 0
    for b in ["pos", "neg"]:
        for seq in test[b]:
            best_match = ""
            best_score = 0
            best_label = ""
            for b2 in ["pos", "neg"]:
                for seq2 in train[b2]:
                    # d -> A dictionary returns the score of any pair of characters.
                    # x -> No gap penalties.
                    al = pairwise2.align.globaldx(seq, seq2, blosum62, score_only = True)
                    print(al)
                    sys.exit(1)

def run_script():
    """Wrapper function"""
    ara_d = unpickle()
    train, test = partition_data(ara_d)
    calc_distance(train, test)

if __name__ == "__main__":
    run_script()