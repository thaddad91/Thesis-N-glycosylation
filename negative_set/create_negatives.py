#!/usr/bin/env python

"""
Author:      Thierry Haddad
Description: Script to remove positive hits from the Arabidopsis proteome
             to create a negative set. Protein-level only.
"""

def run_script():
    # Find all positive labelled proteins
    with open('final_table.txt','r') as f:
        lines = f.readlines()
    pos_ids = [line.split()[0] for line in lines]

    # All proteins from proteome
    with open('ara_proteome.fasta','r') as f:
        lines = f.readlines()
    all_ids = [
            line.split()[0].replace('>','')
            for line in lines if line.startswith('>')
            ]

    # All proteins - positive proteins = negative proteins
    neg_ids = [id for id in all_ids if id not in pos_ids]
    print(len(pos_ids), len(neg_ids), len(all_ids))
    sets = {'pos_ids': pos_ids, 'neg_ids': neg_ids, 'all_ids': all_ids}
    for set in sets:
        with open('{}.txt'.format(set),'w') as f:
            for line in sets[set]:
                f.write(line+'\n')

if __name__ == "__main__":
    run_script()
