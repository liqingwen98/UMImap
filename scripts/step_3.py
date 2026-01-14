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
input_sam = os.path.join(args.path, 'align_transcriptome.sam')
output_sam = os.path.join(args.path, 'align_transcriptome_new.sam')
print(tmp_path)
print(input_sam)
print(output_sam)
for tsv_name in os.listdir(os.path.join(tmp_path, 'tsv')):
    map_umi = {}
    query_umi = read_fasta_with_biopython(os.path.join(tmp_path, 'query', tsv_name.replace('.tsv', '.fasta')))
    ref_umi = read_fasta_with_biopython(os.path.join(tmp_path, 'ref', tsv_name.replace('.tsv', '.fasta')))
    df_umi = read_tsv_with_pandas(os.path.join(tmp_path, 'tsv', tsv_name))
    df_umi = df_umi[df_umi.iloc[:, 2] >= 0]
    dict_umi_ref = dict(zip(df_umi.iloc[:, 0], df_umi.iloc[:, 1]))
    for k, v in dict_umi_ref.items():
        map_umi[query_umi[k]] = ref_umi[v]
    transcript_dict[tsv_name[:-4]] = map_umi

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
