#!/usr/bin/env python3

"""
Author:      Thierry Haddad
Description: Sample 4400 entries (2200 pos/neg) and create sequence-derived
             features tables.
"""

from proteinko import Proteinko
from multiprocessing import Pool

def read_table():
    """Read tab-del table of pos and neg sequons."""
    loc = "/mnt/scratch/hadda003/negative_set/pos_neg_table.txt"
    with open(loc, 'r') as f:
        data = f.readlines()
    create_tables(data)


def create_tables(data = []):
    """
    For all the pos and neg sequons, calculate their sequence-derived
    characteristics like:
    - hydropathy
    - acceptors
    - donors
    - iso-electric point
    - volume
    For every feature, a separate table is created. These will
    afterwards be used in R for ML purposes.
    """
    # Create 4 workers
    pool = Pool(processes = 4)
    # Map calc_dist function to all sequences in the data set
    results = pool.map(calc_dist, data)
    print(len(results))
    with open("table_combined", "w") as f:
        for line in results:
            f.write(line)
    print("Done")

def calc_dist(line):
    """Calculates the distributions per scheme per amino-acid."""
    prt = Proteinko()
    schemes = prt.get_schemes()

    line = line.split("\t")
    seq = line[3].strip()
    if not len(seq) == 15:
        continue
    label = line[1].strip()
    new_line = label
    values = []
    for scheme in schemes:
        # This step will take long for all entries
        dist = list(prt.get_dist(seq, scheme))
        # Check for only floats
        for val in dist:
            try:
                num = float(val)
            except:
                print(seq, val)
                continue  # Skip due to R issues
        dist = list(map(str, dist))  # Writeable
        values.append(dist)

    if len(values) == 6800:
        new_line += "\t" + "\t".join(values) + "\n"
        return new_line
    else:
        print(len(values)

def run_script():
    read_table()

if __name__ == "__main__":
    run_script()
