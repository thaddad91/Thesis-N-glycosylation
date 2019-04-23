#!/usr/bin/env python

"""
Author:     Thierry Haddad
Description: Script to calculate the distance between two A.A. sequences using
             a substitution matrix. Afterwards, assign the class of the lowest
             distance to the new sequence, and calculate accuracy.
"""

from Bio.SubsMat.MatrixInfo import blosum62
from Bio import pairwise2
from sklearn.model_selection import KFold
import pickle
import sys

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

def partition_data(ara_d):
    """Create data partitions for train/test sets."""
    pos_set = []
    neg_set = []
    for id in ara_d:
        for b in ["pos", "neg"]:
            for seq in ara_d[id][b]:
                if seq[1] == "" or len(seq[1]) < 15:
                    continue  # Pass empty sequons
                # Append sequon to positive or negative set
                if b == "pos":
                    pos_set.append(seq[1])
                else:
                    neg_set.append(seq[1])
    #pos_set = pos_set[:100]  # Debugging
    sub_neg = neg_set[:2200]  # Subset for debugging purposes

    kfold = KFold(10, True, 1)


    for train, test in kfold.split(pos_set+sub_neg):
        print(len(train), len(test))
        sys.exit(1)

        # Train (75%) and test (25%) sets
        # train = {
        #     "pos": pos_set[:int(len(pos_set) * 0.9)],
        #     "neg": sub_neg[:int(len(sub_neg) * 0.9)],
        # }
        # test = {
        #     "pos": pos_set[int(len(pos_set) * 0.9):],
        #     "neg": sub_neg[int(len(sub_neg) * 0.9):],
        # }
        #return train, test

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
                    #al = pairwise2.align.globalds(seq, seq2, blosum62, score_only = True)
                    score = 0
                    align = zip(seq, seq2)
                    for pair in align:
                        try:
                            score += blosum62[pair]
                        except:  # Reverse order
                            score += blosum62[(pair[1], pair[0])]
                    if score > best_score:  # Update if higher score
                            best_score = score
                            best_label = b2
                            #best_match = seq2  # Unused so far
            # If best training sequence label equals test label
            if best_label == b:
                good += 1
            else:
                bad += 1
    print(good, bad, good/(good+bad)*100)

def run_script():
    """Wrapper function"""
    ara_d = unpickle()
    train, test = partition_data(ara_d)
    #calc_distance(train, test)

if __name__ == "__main__":
    run_script()