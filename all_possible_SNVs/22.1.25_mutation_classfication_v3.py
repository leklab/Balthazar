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


# make a mutation dictionary

mut_dict={

  "A": {
    "T","G","C"
    },
  "T": {
    "A","G","C"
    },
  "G": {
    "A","C","T"
    },
  "C": {
    "A","G","T"
    }
}

import sys

# input the files
# inFile is the fasta file of the target region
# with Line 0 being the header and Line 1 being the sequence
# The sequence starts from the start codon to the last codon before the stop codon

inFile=open(sys.argv[1],'rt')
outFile=open(sys.argv[2],'wt')

# read the sequence
# sequence[0] is the header
# sequence[1] is the actual sequence

sequence= [line.rstrip('\n') for line in inFile]

# capitalize the sequence
seq=""
for nt in range(0,len(sequence[1])):
    if (sequence[1][nt]=="A") or (sequence[1][nt]=="a"):
        seq=seq+"A"
    if (sequence[1][nt]=="T") or (sequence[1][nt]=="t"):
        seq=seq+"T"
    if (sequence[1][nt]=="G") or (sequence[1][nt]=="g"):
        seq=seq+"G"
    if (sequence[1][nt]=="C") or (sequence[1][nt]=="c"):
        seq=seq+"C"

lenseq=len(seq)

# make a header
outFile.write("site"+'\t'+"codon"+'\t'+"WT_nt"+'\t'+"Variant_nt"+'\t'+"classification"+'\n')


# all start codon mutations are classified as non-start mutations:
codon=1
for ntsite in range(0,3):
    for variant_nt in mut_dict[sequence[1][ntsite]]:
        outFile.write(str(ntsite+1)+'\t'+str(codon)+'\t'+sequence[1][ntsite]+'\t'+variant_nt+'\t'+"Start-loss"+'\n')


# define condons

for codon in range(1,(lenseq/3)):
    codon_seq=seq[codon*3:((codon*3)+3)]

    for codon_nt in range(0,3):


# refresh the mutated codon sequence
        codon_mutated=list(codon_seq)


# generate the mutation
        for x in mut_dict:
            for y in mut_dict[x]:
                if codon_seq[codon_nt]==x:
                    codon_mutated[codon_nt]=y
# classify the mutations

                    if AA["".join(codon_mutated)]=="*":
                        outFile.write(str(codon*3+codon_nt+1)+'\t'+str(codon+1)+'\t'+x+'\t'+y+'\t'+"Nonsense"+'\n')


                    elif AA["".join(codon_mutated)]==AA[codon_seq]:
                        outFile.write(str(codon*3+codon_nt+1)+'\t'+str(codon+1)+'\t'+x+'\t'+y+'\t'+"Synonymous"+'\n')


                    else:
                        outFile.write(str(codon*3+codon_nt+1)+'\t'+str(codon+1)+'\t'+x+'\t'+y+'\t'+"Missense"+'\n')









inFile.close()
outFile.close()
