# Python 3.6 or later
# BioPython 1.78 or later and argparse
import random
from Bio import SeqIO
from Bio.Seq import Seq
import argparse
import csv
import sys

# generate reverse complementary sequence
def rev(seq):
    return str(Seq(seq).reverse_complement())


# make a mutation dictionary
mut_dict = {
    "A": {"T", "G", "C"},
    "T": {"A", "G", "C"},
    "G": {"A", "C", "T"},
    "C": {"A", "G", "T"}
}

degenerate_mut_dict = {
    "A": "B",
    "T": "V",
    "G": "H",
    "C": "D"
}
def read_enzymes_from_file(file_path):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter='\t')
        return {rows['enzyme']: rows for rows in reader}

def main():
    parser = argparse.ArgumentParser(description="Oligo design")
    parser.add_argument("-f", "--file_path", type=str, required=True,
                        help="Path to the input FASTA file.")
    parser.add_argument("-o", "--output", type=str, required=True,
                        help="Path to the output file")
    parser.add_argument("-s", "--start_codon_position", type=int, required=True,
                        help="Position of the first nucleotide of the start codon (counting from 0)")
    parser.add_argument("-e", "--stop_codon_position", type=int, required=True,
                        help=" Position of the first nucleotide of the stop codon (counting from 0)")
    parser.add_argument("-efile", "--enzyme_file", type=str, default="default_enzymes.tsv",
                        help="Path to the enzyme file (TSV)")
    parser.add_argument("-l", "--left", type=int, default=25,
                        help="Length of the upstream sequence of the variant in the oligo")
    parser.add_argument("-r", "--right", type=int, default=19,
                        help="Length of the downstream sequence of the variant in the oligo")
    parser.add_argument("-b", "--barcode_file", type=str, required=True,
                        help="Path to the barcode file ")
    parser.add_argument("-blkl", "--block_length_upper_limit", type=int, default=250,
                        help="Block length upper limit")
    parser.add_argument("-m", "--mode", type=int, choices=[1, 2, 3], required=True,
                        help="Mode of saturation mutagenesis")
# mode1: mutate each nt into the 3 other nts as individual oligos.
# mode2: mutate each nt into the corresponding degenerate nt as 1 oligo.
# mode3: mutate each codon into NNK as 1 oligo.

    args = parser.parse_args()
    try:
        fa_seq = SeqIO.read(args.file_path, 'fasta')
    except Exception as e:
        print("Error reading the file: {}".format(e))
        return

# get gene length
    start_position = args.start_codon_position-args.left
    stop_position = args.stop_codon_position+args.right
    left_gene_seq_right = str(fa_seq.seq[start_position :stop_position]).upper()
    gene_sequence=str(fa_seq.seq[args.start_codon_position :args.stop_codon_position]).upper()
    gene_length = args.stop_codon_position - args.start_codon_position

# decide block number and block length according to the mode
    if args.mode == 3 :
        if gene_length % 3 == 0:
            AAS_length = gene_length // 3
            AAS_block_limit = args.block_length_upper_limit // 3
            AAS_last_block_temp = AAS_length % AAS_block_limit
            block_num_temp = AAS_length // AAS_block_limit
            last_block_temp = AAS_last_block_temp*3

            for i in range(1, AAS_block_limit):
                if i * block_num_temp + AAS_last_block_temp > AAS_block_limit - i:
                    contri = i
                    break

            block_num = block_num_temp + 1
            print("Block number: " + str(block_num))
            block_size_exclude_last_block = (AAS_block_limit - contri + 1) * 3
            print("Block size except the last: " + str(block_size_exclude_last_block))
            block_size_last_block = ((contri - 1) * block_num_temp + AAS_last_block_temp) * 3
            print("Last block size: " + str(block_size_last_block))

        else:
            print("Error: Mode3 is for CDS sequences whose length should be an integer multiple of 3")
            sys.exit(1)

    else:
        last_block_temp = gene_length % args.block_length_upper_limit
        block_num_temp = gene_length // args.block_length_upper_limit
        for i in range(1, args.block_length_upper_limit):
            if i * block_num_temp + last_block_temp > args.block_length_upper_limit - i:
                contri = i
                break

        block_num = block_num_temp + 1
        print("Block number: "+str(block_num))
        block_size_exclude_last_block = args.block_length_upper_limit - contri + 1
        print("Block size except the last: "+str(block_size_exclude_last_block))
        block_size_last_block = (contri - 1) * block_num_temp + last_block_temp
        print("Last block size: "+str(block_size_last_block))

    try:
        enzymes = read_enzymes_from_file(args.enzyme_file)
    except (FileNotFoundError, ValueError) as e:
        print(f"Error reading enzyme file: {e}. Using default enzymes from 'default_enzymes.tsv'.")
        enzymes = read_enzymes_from_file('default_enzymes.tsv')

    for enzyme, enzyme_info in enzymes.items():
        recog_sequence = enzyme_info['recognition_sequence']
        recog_right_size = int(enzyme_info['recog_right_size'])
        insulator_size = int(enzyme_info['insulator_size'])
        if recog_sequence not in left_gene_seq_right and recog_sequence not in rev(left_gene_seq_right):
            type2s_recog = recog_sequence
            print("Type2S enzyme: "+enzyme)
            break
    else:
        print(
            "Try to find another Type IIs enzyme and modify this script or consider introducing synonymous variants to the template")
        type2s_recog = None
        return

    barcode_list = []
    try:
        with open(args.barcode_file, 'r') as file:
            barcode_list = [line.strip() for line in file.readlines()]
    except Exception as e:
        print(f"An error occurred while reading the barcode file: {e}")
        return

    block_adaptor = ['' for _ in range(block_num)]

# get the sequence for each block
# the sequence is with upstream overhang (default 25 nt) and downstream overhang (default 19 nt)
    blocks = [left_gene_seq_right[i-args.left:i + block_size_exclude_last_block+args.right] for i in
              range(args.left, gene_length - block_size_last_block, block_size_exclude_last_block)]
    if last_block_temp > 0:
        blocks.append(left_gene_seq_right[-1 * (block_size_last_block+args.left+args.right):])

# choose the recog_right sequence and make sure the insulator sequence does not present in the corresponding block

    recog_right_memory = []
    recog_counter = 0

    for i in range(block_num):
        while True:

            recog_right = barcode_list[len(barcode_list)-recog_counter-1][0:recog_right_size]
            recog_counter+=1
            check_seq = recog_right[recog_right_size - insulator_size:recog_right_size]

            if (check_seq not in blocks[i] and
                    check_seq not in rev(blocks[i]) and
                    recog_right not in recog_right_memory):
                recog_right_memory.append(recog_right)
                block_adaptor[i] = barcode_list[i] + type2s_recog + recog_right
                break

# Generate mutations
        def mut_1():
            flag = 0
            outFile = open(args.output, 'wt')  # Specify the output file name

            for i in range(block_num):
                for j in range(args.left, len(blocks[i]) - args.right):
                    for x in mut_dict:
                        if blocks[i][j] == x:
                            sorted_values = sorted(mut_dict[x])
                            for y in sorted_values:
                                mutseq = block_adaptor[i] + "".join(rev(blocks[i][j + 1:j + args.right + 1])) + rev(
                                    y) + "".join(rev(blocks[i][j - args.left:j]))
                                flag += 1
                                outFile.write(">" + str(flag) + '\n' + mutseq + '\n')

            outFile.close()

        def mut_2():
            flag = 0
            outFile = open(args.output, 'wt')  # Specify the output file name

            for i in range(block_num):
                for j in range(args.left, len(blocks[i]) - args.right):
                    for x in degenerate_mut_dict:
                        if blocks[i][j] == x:
                            mutseq = (block_adaptor[i] + "".join(rev(blocks[i][j + 1:j + args.right + 1])) +
                                      degenerate_mut_dict[rev(x)] +
                                      "".join(rev(blocks[i][j - args.left:j])))
                            flag += 1
                            outFile.write(">" + str(flag) + '\n' + mutseq + '\n')

            outFile.close()

        def mut_3():
            flag = 0
            outFile = open(args.output, 'wt')  # Specify the output file name

            for i in range(block_num):
                for j in range(args.left, len(blocks[i]) - args.right, 3):
                    if j + 2 < len(blocks[i]) - args.right:
                        mut_block = 'MNN'
                        mutseq = (block_adaptor[i] + "".join(rev(blocks[i][j + 1+2:j + args.right + 1+2])) + mut_block +
                                  "".join(rev(blocks[i][j - args.left:j])))

                        flag += 1
                        outFile.write(">" + str(flag) + '\n' + mutseq + '\n')

            outFile.close()

        if args.mode == 1:
            mut_1()
        elif args.mode == 2:
            mut_2()
        elif args.mode == 3:
            mut_3()

if __name__ == "__main__":
    main()
