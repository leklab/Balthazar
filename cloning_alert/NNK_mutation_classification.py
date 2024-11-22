#make a codon table

AA={}
AA["TTT"]="F"
AA["TTC"]="F"
AA["TTA"]="L"
AA["TTG"]="L"
AA["CTT"]="L"
AA["CTC"]="L"
AA["CTA"]="L"
AA["CTG"]="L"
AA["ATT"]="I"
AA["ATC"]="I"
AA["ATA"]="I"
AA["ATG"]="M"
AA["GTT"]="V"
AA["GTC"]="V"
AA["GTA"]="V"
AA["GTG"]="V"
AA["TCT"]="S"
AA["TCC"]="S"
AA["TCA"]="S"
AA["TCG"]="S"
AA["CCT"]="P"
AA["CCC"]="P"
AA["CCA"]="P"
AA["CCG"]="P"
AA["ACT"]="T"
AA["ACC"]="T"
AA["ACA"]="T"
AA["ACG"]="T"
AA["GCT"]="A"
AA["GCC"]="A"
AA["GCA"]="A"
AA["GCG"]="A"
AA["TAT"]="Y"
AA["TAC"]="Y"
AA["TAA"]="*"
AA["TAG"]="*"
AA["CAT"]="H"
AA["CAC"]="H"
AA["CAA"]="Q"
AA["CAG"]="Q"
AA["AAT"]="N"
AA["AAC"]="N"
AA["AAA"]="K"
AA["AAG"]="K"
AA["GAT"]="D"
AA["GAC"]="D"
AA["GAA"]="E"
AA["GAG"]="E"
AA["TGT"]="C"
AA["TGC"]="C"
AA["TGA"]="*"
AA["TGG"]="W"
AA["CGT"]="R"
AA["CGC"]="R"
AA["CGA"]="R"
AA["CGG"]="R"
AA["AGT"]="S"
AA["AGC"]="S"
AA["AGA"]="R"
AA["AGG"]="R"
AA["GGT"]="G"
AA["GGC"]="G"
AA["GGA"]="G"
AA["GGG"]="G"

import sys

# input the files
# inFile is the fasta file of the target region
# with Line 0 being the header and Line 1 being the sequence
# The sequence starts from the start codon to the last codon before the stop codon

inFile=open(sys.argv[1],'rt')
outFile=open(sys.argv[2],'wt')
start_codon_position = int(sys.argv[3])
stop_codon_position = int(sys.argv[4])

# read the sequence
sequence = ""
next(inFile)
for line in inFile:
    line = line.strip()
    sequence += line
# capitalize the sequence
seq=''
for nt in range(0, len(sequence)):
    if (sequence[nt] == "A") or (sequence[nt] == "a"):
        seq += "A"
    elif (sequence[nt] == "T") or (sequence[nt] == "t"):
        seq += "T"
    elif (sequence[nt] == "G") or (sequence[nt] == "g"):
        seq += "G"
    elif (sequence[nt] == "C") or (sequence[nt] == "c"):
        seq += "C"

lenseq=len(seq[start_codon_position:stop_codon_position])

# make a header
outFile.write("codon" + '\t' + "WT_codon" + '\t' + "Variant_codon" + '\t' + "classification" + '\n')


# all start codon mutations are classified as non-start mutations:
start_codon_seq = seq[start_codon_position:start_codon_position+3]
for first_nt in ['A', 'T', 'G', 'C']:
    for second_nt in ['A', 'T', 'G', 'C']:
        for third_nt in ['G', 'T']:
            variant_codon = first_nt + second_nt + third_nt
            if variant_codon == start_codon_seq:
                outFile.write('1\t' + start_codon_seq + '\t' + variant_codon + '\t' + "Synonymous" + '\n')
            else:
                outFile.write('1\t' + start_codon_seq + '\t' + variant_codon + '\t' + "Start-loss" + '\n')


# define condons
for codon in range(1,(lenseq//3)):
    codon_seq=seq[start_codon_position+codon*3:(start_codon_position+(codon*3)+3)]

# generate and classify the mutation
    for first_nt in ['A', 'T', 'G', 'C']:
        for second_nt in ['A', 'T', 'G', 'C']:
            for third_nt in ['G', 'T']:
                variant_codon = first_nt + second_nt + third_nt
                if AA[variant_codon] == "*":
                    outFile.write(
                        str(codon + 1) + '\t' + codon_seq + '\t' + variant_codon + '\t' + "Nonsense" + '\n')
                elif AA[variant_codon] == AA[codon_seq]:
                    outFile.write(
                        str(codon + 1) + '\t' + codon_seq + '\t' + variant_codon + '\t' + "Synonymous" + '\n')
                else:
                    outFile.write(
                        str(codon + 1) + '\t' + codon_seq + '\t' + variant_codon + '\t' + "Missense" + '\n')

inFile.close()
outFile.close()
