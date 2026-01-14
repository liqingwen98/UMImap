import pysam, os
from collections import defaultdict
from collections import Counter
from Bio import SeqIO
from tqdm import tqdm
import subprocess
from concurrent.futures import ProcessPoolExecutor

def run_command(tran_id, tmp_path):
    cmd = [
        'cd-hit-est',
        '-i', os.path.join(tmp_path, f'{tran_id}'),
        '-o', os.path.join(os.path.dirname(tmp_path), 'cluster', f'{tran_id}'[:-6]),
        '-c', '0.8',
        '-n', '2', 
        '-d', '0', 
        '-M', '100000', 
        '-T', '4',
    ]
    try:
        with open(os.devnull, 'w') as devnull:
            result = subprocess.run(
                cmd,
                stdout=devnull,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        return tran_id, result.stderr
    except subprocess.CalledProcessError as e:
        return tran_id, f"Error: {e.stderr}"

def parallel_command(more_transcript_keys, tmp_path, max_workers=4):
    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(run_command, tran_id, tmp_path) for tran_id in more_transcript_keys]
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
print(tmp_path)
if not os.path.exists(os.path.join(tmp_path, 'cluster')):
    os.makedirs(os.path.join(tmp_path, 'cluster'))
tmp_path = os.path.join(tmp_path, 'ref')
parallel_command(os.listdir(tmp_path), tmp_path, max_workers=args.threads//4)