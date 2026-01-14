from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import os

def remove_umi_from_sequence(sequence, umi, right_flanking_length):
    sequence_str = str(sequence).upper()
    umi = umi.upper()
    start_index = sequence_str.find(umi)
    
    return sequence_str[start_index + len(umi) + right_flanking_length :]
    # if umi in sequence_str:
    #     return sequence_str.replace(umi, "")
    # return sequence_str  

def process_single_umi_fastq(input_fastq, output_fasta, right_flanking):
    updated_records = []

    for record in SeqIO.parse(input_fastq, "fasta"):
        parts = record.id.rsplit("_", 1)  
        if len(parts) != 2:
            print(f"warining: cannot extract UMI from ID {record.id}, skip this record")
            continue

        read_id, umi = parts[0], parts[1]

        updated_seq = remove_umi_from_sequence(record.seq, umi, len(right_flanking))

        new_record = SeqRecord(
            Seq(updated_seq),
            id=read_id+'_'+umi,
            description=read_id+'_'+umi
        )

        updated_records.append(new_record)

    SeqIO.write(updated_records, output_fasta, "fasta")

    print(f"completed!")
    print(f"updated sequences save to: {output_fasta}")
    print(f"total processed {len(updated_records)} sequences")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta", required=True)
    parser.add_argument("-r", "--right_flanking", required=True)
    args = parser.parse_args()

    input_dir = os.path.dirname(args.fasta)

    output_fasta_file = os.path.join(input_dir, "updated_sequences.fasta")

    process_single_umi_fastq(args.fasta, output_fasta_file, args.right_flanking)