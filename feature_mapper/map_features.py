#!/usr/bin/env python3

"""
Author:      Thierry Haddad
Description: Generic script used to map additional tab-delimited feature matrices to the current
             dataset. Might also serve as a template for more complex matrix merges.
"""

import sys
import pickle

# Base data dir
loc = "/mnt/scratch/hadda003/"

def unserialize_ara_d():
    """Test to check serialization integrity."""
    with open(loc+"negative_set/pos_neg_glycosites","rb") as f:
        new_d = pickle.load(f)
    if type(new_d) == dict:
        print("Serialization succes!")
    else:
        print("Serialization failed!")
        print(type(new_d))
        sys.exit(1)
    return new_d



def run_script():
    """Wrapper function"""
    ara_d = unserialize_ara_d()

if __name__ == "__main__":
    run_script()