
TRANSCRIPTOME_FASTA="path_to_reference_transcriptome"
READS_FASTX="path_to_input_read_fastx"
UMImap_DIR="path_to_UMImap_directory"
OUTPUT_DIR="path_to_output_directory"
FORWARD_ADAPTER="sequence_of_forward_adapter"
REVERSE_ADAPTER="sequence_of_reverse_adapter"
LEFT_FLANKING="sequence_of_left_flanking"
RIGHT_FLANKING="sequence_of_right_flanking"
UMI_LENGTH=length_of_UMI
THREADS=number_of_threads

mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"
mkdir -p tmp
cd tmp

python $UMImap_DIR/preprocess/reverse_complement.py --seq1 $FORWARD_ADAPTER --seq2 $REVERSE_ADAPTER --fastx $READS_FASTX --dir $PWD

flexiplex -x "$LEFT_FLANKING"???? -u $(printf '?%.0s' $(seq 1 "$UMI_LENGTH")) -b "" -x ????"$RIGHT_FLANKING" -k "?" -f 7 -e 1 -p $THREADS -n ./flexiumi convert.fasta > flexiumi.fasta

python $UMImap_DIR/preprocess/separate_seq.py --fasta convert.fasta --txt flexiumi_reads_barcodes.txt 

python $UMImap_DIR/preprocess/remove_umi.py --fasta single_umi.fasta --right_flanking $RIGHT_FLANKING

minimap2 -ax map-ont --secondary=no $TRANSCRIPTOME_FASTA updated_sequences.fasta -t $THREADS -u f> align_transcriptome.sam

python $UMImap_DIR/scripts/step_1.py --path $PWD -t $THREADS

python $UMImap_DIR/scripts/step_2.py --path $PWD

python $UMImap_DIR/scripts/step_3.py --path $PWD

python $UMImap_DIR/scripts/step_4.py --path $PWD -t $THREADS

python $UMImap_DIR/scripts/step_5.py --path $PWD

samtools view -bS align_transcriptome_new_cluster.sam > temp.bam

samtools fasta temp.bam > ../final.fasta

cd ..
rm -r tmp

echo "Final UMI detected sequences are in "$PWD"/final.fasta"