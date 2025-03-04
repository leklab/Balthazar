**Saturation Mutagenesis-Reinforced Functional Assays (SMuRF)** is a flexible workflow for deep mutational scanning (DMS) studies. This script is designed to generate mutated oligonucleotide sequences for **PALS-C**.

# Requirements

- Python 3.6 or later
- BioPython 1.78 or later
- argparse
- bwa/0.7.17
- samtools/1.16
- picard/2.27.5


## BioPython Required   

To use this script, ensure you have BioPython installed. You can install BioPython using pip:

    pip install biopython


# Usage

The script is run from the command line and requires several arguments to function. Here is a brief description of each argument:

- -f, --file_path: Path to the input FASTA file.
- -o, --output: Path to the output file.
- -s, --start_codon_position: Position of the first nucleotide of the start codon (counting from 0).
- -e, --stop_codon_position: Position of the first nucleotide of the stop codon (counting from 0).
- -efile, --enzyme_file: Path to the enzyme file (TSV). Defaults to "default_enzymes.tsv".
- -l, --left: Length of the upstream sequence of the variant in the oligo. Defaults to 25.
- -r, --right: Length of the downstream sequence of the variant in the oligo. Defaults to 19.
- -b, --barcode_file: Path to the barcode file.
- -blkl, --block_length_upper_limit: Block length upper limit. Defaults to 250.
- -m, --mode: Mode of saturation mutagenesis. Choices are 1, 2, or 3.

## Modes

- Mode 1: Mutate each nucleotide into the 3 other nucleotides as individual oligos.
- Mode 2: Mutate each nucleotide into the corresponding degenerate nucleotide as 1 oligo.
- Mode 3: Mutate each codon into NNK as 1 oligo.

## Example Command

    python oligo_generator.py -f test_input.fa  -o output.fasta -s 7311 -e 8796 -efile default_enzymes.tsv -l 25 -r 19 -b barcodes.txt -blkl 250 -m 1

# Description

The script performs the following steps:

- Reads the input FASTA file and extracts the sequence.
- Calculates the block number and block length based on the mode and the gene length.
- Reads the enzyme file to find a suitable Type IIs enzyme.
- Reads the barcode file to generate unique adaptors for each oligo.
- Generates mutations based on the chosen mode and writes the results to the output file.

# Notes

## test_input.fa

test_input.fa is the sequence of Addgene plasmid 205150 in the FASTA format.


## output.fasta

output.fasta contains sequences of the oligos to order.

For instance, the example commands generate 4,455 oligos to order.

The synthesized oligos are ready for use for PALS-C cloning.

## Output printed on the screen of the example command:

    Block number: 6
    Block size except the last: 248
    Last block size: 245
    Type2S enzyme: BsmBI
    Block1-specific R1 primer: 5'-CTGACTCACGTCTCCTTAC-3'
    Block2-specific R1 primer: 5'-AGGCGACTCGTCTCGTAGC-3'
    Block3-specific R1 primer: 5'-CGCTTGAACGTCTCCCTAC-3'
    Block4-specific R1 primer: 5'-GCATGAGGCGTCTCGTTCA-3'
    Block5-specific R1 primer: 5'-GTTTCCCGCGTCTCCCTAA-3'
    Block6-specific R1 primer: 5'-TGAAAGACCGTCTCCGATA-3'

## Details of default_enzymes.tsv

We picked 5 Type IIs enzymes and added them in the default_enzymes.tsv file.

|enzyme|recognition_sequence|recog_right_size|insulator_size
|-|-|-|-
|BsmBI|CGTCTC|5|4
|BsaI|GGTCTC|5|4
|BbsI|GAAGAC|6|4
|BspMI|ACCTGC|8|4
|PaqCI|CACCTGC|8|4

The inclusion criteria were:

- Keeping only one enzyme of the isoschizomers.
- Removing enzymes with  a < 4nt stickty-end overhang size.
- Removing enzymes with > 1 cutting sites.
- Removing enzymes with a < 6 nt recognition site length.
- Removing enzymes with a > 8 nt cutting site length.

Other enzymes that meet the inclusion criteria can be added to the default_enzymes.tsv file.

## Barcodes

barcodes.txt were generated with barcode_generator.py.

New barcodes.txt can be generated with the following command:

    python barcode_generator.py

## Visualization of the oligos

Use the following commands to generate the files for iGV visualization:

    module load bwa/0.7.17
    bwa index test_input.fa
    module load samtools/1.16
    samtools faidx test_input.fa
    bwa mem -M -t 8 -R "@RG\tID:test\tPL:test\tSM:test" test_input.fa output.fasta > output.sam 2> output.bwa.sterr.log
    module load picard/2.27.5
    java -Xmx14G -jar $PICARD/picard.jar SortSam INPUT=output.sam OUTPUT= output.bam SORT_ORDER=coordinate
    samtools index output.bam
