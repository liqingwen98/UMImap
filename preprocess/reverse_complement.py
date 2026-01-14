import argparse
from Bio import SeqIO, Align
from tqdm import tqdm
import os
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('--seq1', required=True,)
    parser.add_argument('--seq2', required=True,)
    parser.add_argument('--fastx', required=True)
    parser.add_argument('--dir', required=True)
    args = parser.parse_args()

    seq1 = args.seq1.upper()
    seq2 = args.seq2[::-1].upper()  

    fastx_file = args.fastx
    output_file = os.path.join(args.dir, "convert.fasta")
    suffix = os.path.splitext(fastx_file)[-1]

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'  
    aligner.match_score = 1
    aligner.mismatch_score = -1  
    aligner.gap_score = -1  
    aligner.open_gap_score = -1 
    aligner.extend_gap_score = -1 

    total_count = 0
    reverse_count = 0
    forward_count = 0

    with open(output_file, "w") as out_f:
        for record in tqdm(SeqIO.parse(fastx_file, suffix.replace('.', ''))):
            total_count += 1
            read_seq = str(record.seq)[:100]  

            alignment1 = aligner.align(read_seq, seq1)
            score1 = alignment1.score if alignment1 else -float('inf')

            alignment2 = aligner.align(read_seq, seq2)
            score2 = alignment2.score if alignment2 else -float('inf')

            if score1 < 13 and score2 < 13:
                continue
            elif score1 > score2:
                read_seq  = str(record.seq)[-100:]
                alignment3 = aligner.align(read_seq, str(Seq(seq2).reverse_complement()))
                if alignment3.score < 10:
                    new_record = SeqRecord(
                                    record.seq,
                                    id=record.id,
                                    description=''
                                        )
                    SeqIO.write(new_record, out_f, "fasta")
                else:
                    new_record = SeqRecord(
                                    record.seq[:-100+alignment3[0].coordinates[0, 0]],
                                    id=record.id,
                                    description=''
                                        )
                    SeqIO.write(new_record, out_f, "fasta")
                forward_count += 1
            elif score1 < score2:
                new_record = SeqRecord(
                                record.seq[alignment2[0].coordinates[0, -1]:].reverse_complement(),
                                id=record.id,
                                description=''
                                    )
                SeqIO.write(new_record, out_f, "fasta")
                reverse_count += 1

    print('Total sequences:', total_count)
    print('Reverse sequences:', reverse_count)
    print('Forward sequences:', forward_count)

if __name__ == "__main__":
    main()