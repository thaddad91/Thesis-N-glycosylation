#!/usr/bin/env python

"""
Author:      Thierry Haddad
Description: Script to combine positive and negative datasets in to a
             complete glycosite-level dataset.
"""

import sys
import re
import pickle
from Bio import SeqIO

def create_window(index, sequence):
    """Finds index of sequon in sequence and creates a 15 residue window."""
    sequon = ''
    if index > 8 and len(sequence) - index > 8:  # Complete 15n window
        sequon = sequence[index - 7:index + 8]
    elif index < 8:                              # Motif at very beginning
        sequon = sequence[:index + 8]
    elif len(sequence) - index < 8:              # Motif at very end
        sequon = sequence[index - 7:]
    return sequon

def find_potentials(motif, sequence, ara_d, id, pos):
    """Find potential glycosites in the given sequence."""
    pot_glycos = re.finditer(motif, sequence)  # All potential glycosites
    for pot_glyco in pot_glycos:
        pot_index = pot_glyco.start(0)
        pot_sequon = create_window(pot_index, sequence)
        if pos:  # If from the positive set
            if [pot_index, pot_sequon] not in ara_d[id]["pos"]:  # If not in positive set...
                ara_d[id]["neg"].append([pot_index, pot_sequon])  # Add to negative
        else:  # Negative set
            if id in ara_d:
                if [pot_index, pot_sequon] not in ara_d[id]["pos"]:  # No duplicate pos
                    if [pot_index, pot_sequon] not in ara_d[id]["neg"]:  # No duplicate neg
                        ara_d[id]["neg"].append([pot_index, pot_sequon])
            else:  # New negative
                ara_d[id] = {"pos":[], "sequence":sequence, "neg":[[pot_index, pot_sequon]]}
    return ara_d

def final_table(motif):
    """Positives and negatives from final_table.txt."""
    ara_d = {}
    with open('final_table.txt','r') as f:
        lines = f.readlines()

    not_added = 0
    for line in lines:
        line = line.strip('\n').split('\t')

        # TAIR ID, pre-sequon, protein sequence
        id = line[0]
        pre_seq = line[1].upper()
        sequence = line[2]

        # Clean up sequon
        if '.' in pre_seq:  # Dataset with '.' in sequon from trypsin digestion
            pre_seq = pre_seq.split('.')[1]  # Middle section
        pre_seq = pre_seq.replace('_','')
        motif_loc = re.search(motif, pre_seq).start(0)
        index = sequence.find(pre_seq)  # Find sequon in prot seq
        if index == -1:
            print(id+"\t"+pre_seq+"\t"+str(index))
            not_added += 1
            continue
        index = index + motif_loc  # Find motif in sequon in sequence
        sequon = create_window(index, sequence)  # Get 15n window

        # Add to dict
        if id not in ara_d:
            ara_d[id] = {"pos":[[index, sequon]], "sequence":sequence, "neg":[]}
        else:
            ara_d[id]["pos"].append([index, sequon])

        # Potentials/negatives
        ara_d = find_potentials(motif, sequence, ara_d, id, pos=True)
    print("Not added: ", not_added)
    return ara_d

def neg_from_proteome(motif, ara_d):
    """Negatives from complete Arabidopsis proteome."""
    fasta_sequences = SeqIO.parse(open("ara_proteome.fasta"), "fasta")  # Full proteome
    for fasta in fasta_sequences:
        id, sequence = fasta.id, str(fasta.seq)  # type: (object, str)
        # Find potential glycosites in sequence
        ara_d = find_potentials(motif, sequence, ara_d, id, pos=False)
    return ara_d

def serialize_ara_d(ara_d):
    """Serialize the pos/neg dictionary to a file."""
    with open("pos_neg_glycosites","wb") as f:
        pickle.dump(ara_d, f)
    return True

def unserialize_ara_d():
    """Test to check serialization integrity."""
    with open("pos_neg_glycosites","rb") as f:
        new_d = pickle.load(f)
    if type(new_d) == dict:
        print("Serialization succes!")
    else:
        print("Serialization failed!")
        print(type(new_d))

def final_dict(ara_d):
    """Dataset with similar entries removed and curated."""
    final_d = {}
    c = 0  # Counter for window intersect
    ids = ara_d.keys()
    # Merge all basic TAIR with .1 version
    for id in ids:
        if '.' not in id:
            # Potentials duplicates
            if id+".1" in ids:  # Needs to be merged
                # Same sequence?
                if ara_d[id]["sequence"] == ara_d[id+".1"]["sequence"]:
                    # Same sequons?
                    set1_p = ara_d[id]["pos"]
                    set1_n = ara_d[id]["neg"]
                    set2_p = ara_d[id+".1"]["pos"]
                    set2_n = ara_d[id+".1"]["neg"]
                    set_p = set1_p+set2_p
                    set_n = set1_n+set2_n

                    # Remove neg duplicates + window intersect
                    n2, p, c = remove_neg_dups(set_n, set_p, c)
                    # Update dict
                    final_d[id+".1"] = {"pos":p, "neg":n2, "sequence":ara_d[id]["sequence"]}
            else:  # Doesn't need to be merged
                n2, p, c = remove_neg_dups(ara_d[id]["neg"], ara_d[id]["pos"], c)
                final_d[id] = {"pos": p, "neg": n2, "sequence": ara_d[id]["sequence"]}

    # Now add the rest
    for id in ids:
        if '.' in id:
            if id not in final_d:  # Only new entries
                n2, p, c = remove_neg_dups(ara_d[id]["neg"], ara_d[id]["pos"], c)
                final_d[id] = {"pos": p, "neg": n2, "sequence": ara_d[id]["sequence"]}
    print(c)
    return final_d

def remove_neg_dups(set_n, set_p, c):
    """Function to remove negative duplicates and in-positive-window negatives."""
    p_set = set(map(tuple, set_p))
    p = list(map(list, p_set))
    n_set = set(map(tuple, set_n))
    n = list(map(list, n_set))
    if len(set_p) > 0 and len(set_p) != len(p):
        print(set_p)
        print(p)
        print("\n")

    # Removes negatives that are found in positives
    n2 = n
    for e in n:
        if e in p:  # If negative found in positive
            n2 = list(filter((e).__ne__, n2))
        else:
            # Remove negs from same window as pos
            for pair in p:
                if abs(e[0]-pair[0]) <= 7:
                    #print(e, pair)
                    c += 1
                    n2 = list(filter((e).__ne__, n2))
    return n2, p, c

def create_table(ara_d):
    """Create tab-delimited file of the dataset."""
    with open("pos_neg_table.txt","w") as f:
        for id in ara_d:
            sequence = ara_d[id]["sequence"]
            posneg = ["pos", "neg"]
            for pn in posneg:
                for pair in ara_d[id][pn]:
                    index = pair[0]
                    sequon = pair[1]
                    f.write("{}\t{}\t{}\t{}\t{}\n".format(id, pn, index, sequon, sequence))
    return True

def run_script():
    """Run the entire script."""
    # Nx[S/T] glycosite motif
    motif = re.compile('N[^P][ST]')

    # Final dict with pos/neg per TAIR ID
    ara_d = final_table(motif)

    # Complete proteome that results in negatives
    ara_d = neg_from_proteome(motif, ara_d)
    ara_d = final_dict(ara_d)
    print(len(ara_d))
    pos = sum([len(ara_d[id]["pos"]) for id in ara_d])
    print(pos)
    neg = sum([len(ara_d[id]["neg"]) for id in ara_d])
    print(neg)
    serialize_ara_d(ara_d)
    #unserialize_ara_d()  # Test, not needed normally

    create_table(ara_d)

if __name__ == "__main__":
    run_script()
