#!/usr/bin/env python3

"""
Author:      Thierry Haddad
Description: Sample 4400 entries (2200 pos/neg) and create sequence-derived
             features tables.
"""

from proteinko import Proteinko

prot = Proteinko()

# Add read table func, see negative dir

def create_tables(data = []):
    """
    For all the pos and neg sequons, calculate their sequence-derived
    characteristics like hydropathy, acceptors, donors, iso-electric point
    and volume. For every feature, a separate table is created. These will
    afterwards be used in R for ML purposes.
    """
    schemes = prt.get_schemas()
    # Run for every feature
    for scheme in schemes:
        table = []  # New table
        for line in data:
            seq = data[3]
            label = data[1]
            # This step will take long for all entries
            dist = list(prt.get_dist(seq, scheme))
            new_line = label + "\t" + "\t".join(dist) + "\n"
            table.append(new_line)
        # Write feature table to tab-del file
        if table != []:
            with open("table_"+scheme, "w") as f:
                for line in table:
                    f.write(line)

def run_script():
    pass

if __name__ == "__main__":
    run_script()
