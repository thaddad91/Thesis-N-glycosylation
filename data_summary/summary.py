#!/usr/bin/env python3

# Author: Thierry Haddad
# Description: Parses 4 csv files and gives simple statistics.

import os
import sys
import re

# 4 Main datafile saved as tab-delimited
data_dir = '/home/hadda003/thesis/start_data/tab_del/'
data = []
files = os.listdir(data_dir)
tl_files = [file for file in files if file.endswith('.txt')]

for file in tl_files:
    print(file)

    # Specific identifier and sequence index per file
    if "Ziel" in file:
        id = 0
        seq = 6
    elif "Song" in file:
        id = 0
        seq = 5
    elif "Xu" in file:
        id = 0
        seq = 4
    elif "Zeng" in file:
        id = 0
        seq = 3
    else:
        print("Not recognized file: ", file)
        sys.exit(1)
    
    try:
        with open(data_dir+file, 'r') as df:
            file_data = df.readlines()
    except IOError as e:
        print("Can't open file")
        print(e)

    last_id = ''
    for line in file_data[1:]:
        line = line.split('\t')
        prot_id = line[id].upper()
        # Some entries list multiple sequences per identifier
        if prot_id == '':
            prot_id = last_id
        else:
            last_id = prot_id
        sequence = line[seq]
        data.append([prot_id, sequence])
print("Total combined entries, inc. redundant ones: ", len(data))

# Create non-redundant version
nonrd_data = []
for pair in data:
    if pair not in nonrd_data:
        nonrd_data.append(pair)
nonrd_data.sort()
print("Total combined entries after redundancy check: ", len(nonrd_data))
with open('combined_data.txt','w') as cb:
    for entry in nonrd_data:
        #if 'AT' in entry[0] and '.' in entry[0]:
            #print(entry[0])
        line = entry[0]+"\t"+entry[1]+"\n"
        cb.write(line)
#tair_ids = list(set([entry[0] for entry in nonrd_data]))
#print(tair_ids)
#with open('all_tairs.txt','w') as at:
#    for entry in tair_ids:
#        at.write(entry+'\n')

# Parse representative list from TAIR bulk
repr_dict = {}
with open('final_representatives.txt','r') as rep:
    repr_data = rep.readlines()
for item in repr_data:
    item = item.split(' ')
    repr_dict[item[0]] = item[-1].rstrip('\n')

# Final table of ID+sequon+prot sequence
final_table = []
not_found = 0
for entry in nonrd_data:
    # Regex to find motif
    reg = re.search('N[^P][ST]', str(entry[1].upper()))
    if not reg:
        print("No motif: ", entry[1])
        not_found += 1
        continue
#    try:
#        motif_hits = re.findall('N[^P][ST]', repr_dict[entry[0]])
#    except KeyError as e:
        #print(e)
#        not_found += 1
        #with open('not_found.txt','a') as nf:
            #nf.write(entry[0]+'\n')
#        motif_hits = []
#        pass
    # If any motifs found
    #if len(motif_hits) > 0 :
        #print("Hits: ", entry[0], len(motif_hits))
    if entry[0].upper() in repr_dict.keys():
        line = entry[0]+"\t"+entry[1]+"\t"+repr_dict[entry[0]]+"\n"
        final_table.append(line)
    elif entry[0].upper()+".1" in repr_dict.keys():
        line = entry[0]+"\t"+entry[1]+"\t"+repr_dict[entry[0]+".1"]+"\n"
        final_table.append(line)
    else:
        print("Not found: ", entry[0])
print(not_found)

with open("final_table.txt","w") as ft:
    for line in final_table:
        ft.write(line)

#ids = []
#for entry in nonrd_data:
#    if 'AT' in entry[0] and '-' not in entry[0]:
#        if entry[0] not in ids:
#            ids.append(entry[0])
#
#with open('identifiers.txt','w') as wf:
#    for id in ids:
#        wf.write(id+'\n')
