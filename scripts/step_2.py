#!/usr/bin/env python
import sys, os
from tqdm import tqdm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--path', required=True)
args = parser.parse_args()
tmp_path = os.path.join(args.path, 'tmp')
print(tmp_path)
if not os.path.exists(os.path.join(tmp_path, 'tsv')):
    os.makedirs(os.path.join(tmp_path, 'tsv'))
for input_file_name in tqdm(os.listdir(os.path.join(tmp_path, 'alignments'))):
    input_file = os.path.join(tmp_path, 'alignments', input_file_name)
    output_file = os.path.join(tmp_path, 'tsv', input_file_name.replace('.aln', '.tsv'))

    results = []
    query_id, target_id, identity = None, None, None

    flag = 0
    with open(input_file, 'r') as f:
        if len(f.readlines()) < 4:
            continue
    with open(input_file, 'r') as f:
        for line in tqdm(f):
            line = line.strip()
            if line.startswith("Query >"):
                query_id = line.split('>')[1].strip()
            elif line.startswith("%Id"):
                flag = 1
            elif flag == 1:
                info = line.strip().split()
                target_id = info[-1]
                flag = 0
            elif '|' in line:
                identity = line.count("|")
            elif "cols" in line and query_id and target_id:
                results.append((query_id, target_id, identity))
                query_id, target_id, identity = None, None, None

    with open(output_file, 'w') as out:
        out.write("Query_ID\tTarget_ID\tIdentity(%)\n")
        for record in results:
            out.write(f"{record[0]}\t{record[1]}\t{record[2]}\n")
