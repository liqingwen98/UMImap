import pysam, os
from collections import defaultdict
from collections import Counter
from Bio import SeqIO
from tqdm import tqdm
import subprocess
from concurrent.futures import ProcessPoolExecutor

def extract_umi(read_name):

    parts = read_name.rsplit("_", 1)  
    if len(parts) == 2:
        return parts[1]  
    else:
        return None  

def count_unique_transcripts(sam_file, tmp_path):

    umi_to_transcripts = defaultdict(set)
    transcripts_to_umi = defaultdict(set)
    depth_count = defaultdict(int)
    total = 0

    with pysam.AlignmentFile(sam_file, "r") as sam:
        for read in tqdm(sam):
            if not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                read_name = read.query_name
                transcript_id = read.reference_name  

                umi = extract_umi(read_name)
                if umi is None:
                    print(f"warning: cannot extract UMI from read {read_name}, skip")
                    continue

                umi_to_transcripts[umi].add(transcript_id)
                transcripts_to_umi[transcript_id].add(umi)
                depth_count[umi] += 1
                
                total = total+1

    results = {}
    for umi, transcripts in umi_to_transcripts.items():
        unique_count = len(transcripts)  
        results[umi] = unique_count

    only_one = []
    more_two = []
    for umi, num in results.items():
        if num == 1:
            only_one.append(umi)
        else:
            more_two.append(umi)
            
    only_one_set = set(only_one)
    more_two_set = set(more_two)  
    
    one_transcript = {}
    more_transcript = {}
    for transcripts, umi in tqdm(transcripts_to_umi.items()):
        one_transcript[transcripts] = umi & only_one_set
        more_transcript[transcripts] = umi & more_two_set
    
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)
    if not os.path.exists(os.path.join(tmp_path, 'ref')):
        os.makedirs(os.path.join(tmp_path, 'ref'))
    if not os.path.exists(os.path.join(tmp_path, 'query')):
        os.makedirs(os.path.join(tmp_path, 'query'))
    if not os.path.exists(os.path.join(tmp_path, 'alignments')):
        os.makedirs(os.path.join(tmp_path, 'alignments'))
        
    for tran_id, ref_umis in one_transcript.items():
        with open(os.path.join(tmp_path, 'ref', tran_id+'.fasta'), "w") as f:
            for i, seq in enumerate(ref_umis, start=1):
                f.write(f">ref_{i}\n{seq}\n")
                
    for tran_id, query_umis in more_transcript.items():
        with open(os.path.join(tmp_path, 'query', tran_id+'.fasta'), "w") as f:
            for i, seq in enumerate(query_umis, start=1):
                f.write(f">query_{i}\n{seq}\n")
    
    return one_transcript, more_transcript

def run_vsearch(tran_id, tmp_path):
    cmd = [
        'vsearch',
        '--usearch_global', os.path.join(tmp_path, 'query', f'{tran_id}.fasta'),
        '--db', os.path.join(tmp_path, 'ref', f'{tran_id}.fasta'),
        '--id', '0.0',
        '--alnout', os.path.join(tmp_path, 'alignments', f'{tran_id}.aln'),
        '--maxaccepts', '1',
        '--query_cov', '0.0',
        '--minseqlength', '16',
        '--threads', '4',  
        '--wordlength', '3'
    ]
    try:
        result = subprocess.run(
            cmd,
            stderr=subprocess.PIPE,
            text=True,
            check=True
        )
        return tran_id, result.stderr
    except subprocess.CalledProcessError as e:
        return tran_id, f"Error: {e.stderr}"

def parallel_vsearch(more_transcript_keys, tmp_path, max_workers=4):
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        # Map the run_vsearch function to all transcript IDs
        futures = [executor.submit(run_vsearch, tran_id, tmp_path) for tran_id in more_transcript_keys]
        # Use tqdm to show progress
        for future in tqdm(futures, total=len(more_transcript_keys), desc="Processing transcripts"):
            tran_id, stderr = future.result()
            if "Error" in stderr:
                print(f"Failed for {tran_id}: {stderr}")

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--path', required=True)
parser.add_argument('-t', '--threads', required=False, default=8, type=int)
args = parser.parse_args()
tmp_path =  os.path.join(args.path, 'tmp')
if not os.path.exists(tmp_path):
    os.makedirs(tmp_path)
one_transcript, more_transcript = count_unique_transcripts(os.path.join(args.path, 'align_transcriptome.sam'),
                                                           tmp_path)
parallel_vsearch(more_transcript.keys(), tmp_path, max_workers=args.threads//4)