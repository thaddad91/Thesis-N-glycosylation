#!/usr/bin/env python3

"""
Author:      Thierry Haddad
Description: Sample 4400 entries (2200 pos/neg) and create sequence-derived
             features tables.
"""

from proteinko import Proteinko
from multiprocessing import Process

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
    prt = Proteinko()
    schemes = prt.get_schemas()
    procs = []
    # Run multiprocess for every feature
    for scheme in schemes:
        proc = Process(target=calc_dist, args=(data, scheme, prt))
        procs.append(proc)
        proc.start()

    for proc in procs:
        proc.join()

def calc_dist(data, scheme, prt):
    """Calculates the scheme disytribution per amino-acid."""
    table = []  # New table
    total = len(data)

    for c, line in enumerate(data):
        if not c%1000:
            print(scheme, c, total)  # Debug counter
        line = line.split("\t")
        seq = line[3]
        if not len(seq) == 15:
            continue
        label = line[1]

        # This step will take long for all entries
        dist = list(prt.get_dist(seq, scheme))
        dist = list(map(str, dist))  # Writeable
        new_line = label + "\t" + "\t".join(dist) + "\n"
        table.append(new_line)

    # Write feature table to tab-del file
    if table != []:
        print("Writing table...")
        with open("table_"+scheme, "w") as f:
            for line in table:
                f.write(line)
    print("Done.")

def run_script():
    read_table()

if __name__ == "__main__":
    run_script()
