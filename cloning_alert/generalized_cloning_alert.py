# Python 3.6 or later
# BioPython 1.78 or later
#group1 will introduce a cut before the variants which will remove the variant from the top strand
#group2 will introduce a cut right after the variants which might affect the elongation of the top strand

import sys
from Bio.Seq import Seq
import argparse

def rev(seq):
    return str(Seq(seq).reverse_complement())

def search_rec_sites_SNV(seq, lenseq , rec_sites , group, rec_site_length , inclasslines , outFile):
    for i in range(lenseq - rec_site_length + 1):
        for n in range(rec_site_length):
            if seq[i + n] != rec_sites[n] and seq[i:i + n] + seq[i + n + 1:i + rec_site_length] == rec_sites[0:n] + rec_sites[n + 1:rec_site_length]:
                for j in range(1, len(inclasslines)):
                    inclassthisline = inclasslines[j].split("\t")
                    if i + n == int(inclassthisline[0]) - 1 and rec_sites[n] == inclassthisline[3]:
                        outFile.write(inclasslines[j] + '\t' + group + '\n')

def search_rec_sites_NNK(seq, lenseq , rec_sites , group, rec_site_length , inclasslines , outFile):
   for i in range(lenseq-rec_site_length+1):
        for n in range(rec_site_length):
            if (i + n)%3==0 and seq[i + n:i + n + 3] != rec_sites[n:n + 3] and seq[i:i + n] + seq[i + n + 3:i + rec_site_length] == rec_sites[0:n] + rec_sites[n + 3:rec_site_length]:
                for j in range(1, len(inclasslines)):
                    inclassthisline = inclasslines[j].split("\t")
                    if (i + n)//3 == int(inclassthisline[0])-1 and rec_sites[n:n + 3] == inclassthisline[2]:
                        outFile.write(inclasslines[j] + '\t' + group + '\n')

# input the files
# inFile is the fasta file of the target region
# with Line 0 being the header and Line 1 being the sequence
# The sequence starts from the start codon to the last codon before the stop codon

def main():
    parser = argparse.ArgumentParser(description="Cloning alert")
    parser.add_argument("-f", "--file_path", type=str, required=True,
                        help="Path to the input FASTA file.")
    parser.add_argument("-i", "--inClass", type=str, required=True,
                        help="Path to the inClass file (TXT).")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to the output file (TXT).")
    parser.add_argument("-s", "--start_position", type=int, required=True,
                        help="Position of the first nucleotide of the start codon (counting from 0).")
    parser.add_argument("-e", "--stop_position", type=int, required=True,
                        help="Position of the first nucleotide of the stop codon (counting from 0).")
    parser.add_argument("-efile", "--enzyme_file", type=str, default="default_enzymes.tsv",
                        help="Path to the enzyme file (TSV).")
    parser.add_argument("-n", "--type2s_enzyme_name", type=str, required=True,
                        help="Choice of type2s enzyme.")
    parser.add_argument("-m", "--mode", type=int, choices=[1, 2], required=True,
                        help="Mode of alert check")
    args = parser.parse_args()

    enzyme_settings = {}

    try:
        with open(args.enzyme_file, 'rt') as file:
            next(file)
            for line in file:
                parts = line.strip().split('\t')
                enzyme_name, recognition_sequence, _, _= parts
                enzyme_settings[enzyme_name] = {
                    'sense_rec_site': recognition_sequence,
                    'antisense_rec_site': rev(recognition_sequence),
                    'rec_site_lengths':len(recognition_sequence)
                }
    except (FileNotFoundError, ValueError) as e:
        print(f"Error reading enzyme file: {e}. Using default enzymes from 'default_enzymes.tsv'.")
        try:
            with open('default_enzymes.tsv', 'rt') as file:
                next(file)
                for line in file:
                    parts = line.strip().split('\t')
                    enzyme_name, recognition_sequence, _, _= parts
                    enzyme_settings[enzyme_name] = {
                        'sense_rec_site': recognition_sequence,
                        'antisense_rec_site': rev(recognition_sequence),
                        'rec_site_lengths': len(recognition_sequence)
                    }
        except (FileNotFoundError, ValueError) as e:
            print(f"Error reading default enzyme file: {e}. Exiting.")
            sys.exit(1)

    if args.type2s_enzyme_name in enzyme_settings:
        settings = enzyme_settings[args.type2s_enzyme_name]
        sense_rec_site = settings["sense_rec_site"]
        antisense_rec_site = settings["antisense_rec_site"]
        rec_site_length = settings["rec_site_lengths"]
    else:
        print(f"Try to find another Type IIs enzyme")
        sys.exit(1)


    with open(args.inClass, 'rt') as inClass:
        inclasslines = [line.rstrip('\n') for line in inClass]
# read the sequence
# capitalize the sequence
    with open(args.file_path, 'rt') as inFile:
        sequence = ""
        next(inFile)
        for line in inFile:
            line = line.strip()
            sequence += line
    seq = ''
    for nt in range(0, len(sequence)):
        if (sequence[nt] == "A") or (sequence[nt] == "a"):
            seq += "A"
        elif (sequence[nt] == "T") or (sequence[nt] == "t"):
            seq += "T"
        elif (sequence[nt] == "G") or (sequence[nt] == "g"):
            seq += "G"
        elif (sequence[nt] == "C") or (sequence[nt] == "c"):
            seq += "C"
    seq = seq[args.start_position:args.stop_position]
    lenseq = len(seq)
    if args.mode == 1:
        with open(args.output, 'wt') as outFile:
            outFile.write("site"+'\t'+"codon"+'\t'+"WT_nt"+'\t'+"Variant_nt"+'\t'+"classification"+'\t'+"alert_group"+'\n')
            search_rec_sites_SNV(seq, lenseq , sense_rec_site, "group_2", rec_site_length , inclasslines , outFile)
            search_rec_sites_SNV(seq, lenseq , antisense_rec_site, "group_1", rec_site_length , inclasslines , outFile)

    else:
        with open(args.output, 'wt') as outFile:
            outFile.write("codon"+'\t'+"WT_codon"+'\t'+"Variant_codon"+'\t'+"classification"+'\t'+"alert_group"+'\n')
            search_rec_sites_NNK(seq, lenseq, sense_rec_site, "group_2", rec_site_length, inclasslines, outFile)
            search_rec_sites_NNK(seq, lenseq, antisense_rec_site, "group_1", rec_site_length, inclasslines, outFile)

if __name__ == "__main__":
    main()