import argparse
from Bio import SeqIO
from collections import defaultdict
import os
from tqdm import tqdm

def find_umis_in_sequence(sequence, umi_list):
    sequence = str(sequence).upper()
    found_umis = []
    for umi in umi_list:
        if umi in sequence:
            found_umis.append(umi)
    return found_umis

def process_fastq(input_fastq, input_txt, output_single, output_multiple):
    umi_dict = defaultdict(list) 
    with open(input_txt, 'r') as txt_file:
        next(txt_file)  
        for line in tqdm(txt_file):
            fields = line.strip().split('\t')
            read_id, _, _, _, umi = fields
            umi_dict[read_id].append(umi)

    single_umi_records = []
    multiple_umi_records = []

    for record in tqdm(SeqIO.parse(input_fastq, "fasta")):
        read_id = record.id
        if read_id in umi_dict:
            expected_umis = umi_dict[read_id]  
            found_umis = find_umis_in_sequence(record.seq, expected_umis)

            if len(found_umis) == 1:
                new_id = f"{read_id}_{found_umis[0]}"
                record.id = new_id
                record.name = new_id
                record.description = new_id
                single_umi_records.append(record)
            elif len(found_umis) > 0:
                umi_string = "_".join(found_umis)
                new_id = f"{read_id}_{umi_string}"
                record.id = new_id
                record.name = new_id
                record.description = new_id
                multiple_umi_records.append(record)

    SeqIO.write(single_umi_records, output_single, "fasta")
    SeqIO.write(multiple_umi_records, output_multiple, "fasta")

    print(f"completed!")
    print(f"single umi containing sequences save to: {output_single}")
    print(f"multiple umi containing sequences save to: {output_multiple}")
    print(f"single umi file contains {len(single_umi_records)} sequences")
    print(f"multiple umi file contains {len(multiple_umi_records)} sequences")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-t", "--txt", required=True)
    args = parser.parse_args()

    input_dir = os.path.dirname(args.fasta)

    output_single_umi = os.path.join(input_dir, "single_umi.fasta")
    output_multiple_umi = os.path.join(input_dir, "multi_umi.fasta")

    process_fastq(args.fasta, args.txt, output_single_umi, output_multiple_umi)