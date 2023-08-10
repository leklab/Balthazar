clip1=7
clip2=89


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

in_apSNVs=open(sys.argv[1],'rt')
in_clinvar=open(sys.argv[2],'rt')
in_seq=open(sys.argv[3],'rt')
outFile=open(sys.argv[4],'wt')


# read the sequence


in_apSNVs_lines=[line.rstrip('\n') for line in in_apSNVs]


in_clinvar_lines=[line.rstrip('\n') for line in in_clinvar]
in_seq_lines=[line.rstrip('\n') for line in in_seq]

seq=in_seq_lines[1]

# make a header
outFile.write("site"+'\t'+"codon"+'\t'+"WT_nt"+'\t'+"Variant_nt"+'\t'+"WT_AA"+'\t'+"Variant_AA"+'\t'+"classification"+'\t'+"Clinvar_clinical_significance"+'\n')





for i in range(1, len(in_apSNVs_lines)):
    apSNVs_this_line=in_apSNVs_lines[i].split("\t")
    site=int(apSNVs_this_line[0])
    codon=apSNVs_this_line[1]
    WT_nt=apSNVs_this_line[2]
    Variant_nt=apSNVs_this_line[3]
    classification=apSNVs_this_line[4]
    if site%3==1:
        WT_AA=AA[seq[site-1:site+2]]
        Variant_AA=AA[Variant_nt+seq[site:site+2]]

    if site%3==2:
        WT_AA=AA[seq[site-2:site+1]]
        Variant_AA=AA[seq[site-2]+Variant_nt+seq[site]]

    if site%3==0:
        WT_AA=AA[seq[site-3:site]]
        Variant_AA=AA[seq[site-3:site-1]+Variant_nt]

    clinvarnote="NA"
    for LLL in range(clip1,len(in_clinvar_lines)-clip2):
        cvtl=in_clinvar_lines[LLL].split("\t")
        if len(cvtl)>=5:
            cvnote=cvtl[4].split("(")
            cvtldot=cvtl[0].split("c.")
            if len(cvtldot)>=2:
                cvtlspw=cvtldot[1].split(">")
                if len(cvtlspw)>=2:
                    cvtlvnt=cvtlspw[1].split(" (")
                    if cvtlspw[0]==str(site)+WT_nt and cvtlvnt[0]==Variant_nt:
                        clinvarnote=cvnote[0]
    outFile.write(str(site)+'\t'+str(codon)+'\t'+WT_nt+'\t'+Variant_nt+'\t'+WT_AA+'\t'+Variant_AA+'\t'+classification+'\t'+clinvarnote+'\n')





in_apSNVs.close()
in_clinvar.close()
in_seq.close()
outFile.close()
