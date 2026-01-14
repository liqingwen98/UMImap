from Bio import SeqIO
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys, os, pysam
from tqdm import tqdm

def read_tsv_with_pandas(file_path):
    df = pd.read_csv(file_path, sep='\t', encoding='utf-8')
    return df

def read_fasta_with_biopython(file_path):
    fasta_dict = {}
    for record in SeqIO.parse(file_path, "fasta"):
        fasta_dict[record.id] = str(record.seq)
    return fasta_dict

def extract_umi(read_name):
    parts = read_name.rsplit("_", 1)  
    if len(parts) == 2:
        return parts[1]  
    else:
        return None  

transcript_dict = {}
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--path', required=True)
args = parser.parse_args()
tmp_path = os.path.join(args.path, 'tmp')
input_sam = os.path.join(args.path, 'align_transcriptome_new.sam')
output_sam = os.path.join(args.path, 'align_transcriptome_new_cluster.sam')
print(tmp_path)
print(input_sam)
print(output_sam)
for file_name in tqdm(os.listdir(os.path.join(tmp_path, 'ref'))):
    map_umi = {}
    ref_umi = read_fasta_with_biopython(os.path.join(tmp_path, 'ref', file_name))
    with open(os.path.join(tmp_path, 'cluster', file_name.replace('.fasta', '.clstr')), 'r') as file:
        clusters_data = file.read()
    clusters = clusters_data.strip().split('>Cluster ')[1:]
    
    cluster_dict = {}
    for cluster in clusters:
        lines = cluster.strip().split('\n')
        cluster_num = lines[0].split()[0]
        ref_seq = None
        other_seqs = []
        
        for line in lines[1:]: 
            if '*' in line:
                ref_seq = line.split('>')[1].split('...')[0]  
            else:
                other_seqs.append(line.split('>')[1].split('...')[0])  

        if len(other_seqs) > 0:
            for seq in other_seqs:
                cluster_dict[seq] = ref_seq

    for k, v in cluster_dict.items():
        map_umi[ref_umi[k]] = ref_umi[v]
    transcript_dict[file_name[:-6]] = map_umi

with pysam.AlignmentFile(input_sam, "r") as infile, pysam.AlignmentFile(output_sam, "w", template=infile) as outfile:
    for read in tqdm(infile):
        transcript_id = read.reference_name  
        read_name = read.query_name
        umi = extract_umi(read_name)
        
        if transcript_id in transcript_dict:
            umi_dict = transcript_dict[transcript_id]
            if umi in umi_dict:
                last_underscore = read_name.rfind('_')
                if last_underscore != -1:
                    base_name = read_name[:last_underscore]
                    read.query_name = f"{base_name}_{transcript_dict[transcript_id][umi]}"
        
        outfile.write(read)

print('Finished!')
