This script identifies groups of variants that might be underrepresented in the pool due to the limitation of **PALS-C**.

_group1_: will introduce a cut before the variants which will remove the variant from the top strand  
_group2_: will introduce a cut right after the variants which might affect the elongation of the top strand

# Requirements

- Python 3.6 or later
- FASTA file with DNA sequence
- Classification file with variant information
- BioPython 1.78 or later
- argparse

## BioPython Required   

To use this script, ensure you have BioPython installed. You can install BioPython using pip:

    pip install biopython


# Usage

The script is run from the command line and requires several arguments to function. Here is a brief description of each argument:

- -f, --file_path: Path to the input FASTA file.
- -i, --inClass: Path to the inClass file (TXT).
- -o, --output: Path to the output file.
- -s, --start_position: Position of the first nucleotide of the start codon (counting from 0).
- -e, --stop_position: Position of the first nucleotide of the stop codon (counting from 0).
- -efile, --enzyme_file: Path to the enzyme file (TSV). Defaults to "default_enzymes.tsv".
- -n, --type2s_enzyme_name: Choice of type2s enzyme.
- -m, --mode: Mode of alert check.

## Modes

- Mode 1: Mutate each nucleotide into the 3 other nucleotides as individual oligos.
- Mode 2: Mutate each codon into NNK as 1 oligo.

## Example Command

    python generalized_cloning_alert.py -f test_input.fa -i all_possible_SNVs_test.txt  -o output_test_mode1.txt -s 7311 -e 8796 -efile default_enzymes.tsv -n BsmBI -m 1

    python generalized_cloning_alert.py -f test_input.fa -i all_possible_NNKmut_test.txt  -o output_test_mode2.txt -s 7311 -e 8796 -efile default_enzymes.tsv -n BsmBI -m 2

# Description

The script performs the following steps:

- Reads all the inputs ,and convert DNA sequence to uppercase.
- Search for sense sites and antisense sites.
- Writes the results to the output file.


# Notes

## test.fa

test.fa is the sequence of Addgene plasmid 205150.

## inClassfile.txt

all_possible_SNVs_test.txt were generated with SNV_mutation_classification.py.
all_possible_NNKmut_test.txt were generated with NNK_mutation_classification.py.

### Example Command

    python SNV_mutation_classification.py test_input.fa all_possible_SNVs_test.txt 7311 8796
    python NNK_mutation_classification.py test_input.fa all_possible_NNKmut_test.txt 7311 8796    


## Details of default_enzymes.tsv

We picked 5 Type IIs enzymes and added them in the default_enzymes.tsv file.

The inclusion criteria were:

- Keeping only one enzyme of the isoschizomers.
- Removing enzymes with  a < 4nt stickty-end overhang size.
- Removing enzymes with > 1 cutting sites.
- Removing enzymes with a < 6 nt recognition site length.
- Removing enzymes with a > 8 nt cutting site length.

Other enzymes that meet the inclusion criteria can be added to the default_enzymes.tsv file.

|enzyme|recognition_sequence|recog_right_size|insulator_size
|-|-|-|-
|BsmBI|CGTCTC|5|4
|BsaI|GGTCTC|5|4
|BbsI|GAAGAC|6|4
|BspMI|ACCTGC|8|4
|PaqCI|CACCTGC|8|4
