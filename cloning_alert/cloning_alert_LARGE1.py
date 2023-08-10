
#group1 will introduce a cut before the variants which will remove the variant from the top strand
#group2 will introduce a cut right after the variants which might affect the elongation of the top strand

#type2S enzyme

sense_rec_site="GGTCTC"
antisense_rec_site="GAGACC"


import sys

# input the files
# inFile is the fasta file of the target region
# with Line 0 being the header and Line 1 being the sequence
# The sequence starts from the start codon to the last codon before the stop codon

inFile=open(sys.argv[1],'rt')
inClass=open(sys.argv[2],'rt')
outFile=open(sys.argv[3],'wt')

# read the sequence
# sequence[0] is the header
# sequence[1] is the actual sequence

sequence = [line.rstrip('\n') for line in inFile]
inclasslines =[line.rstrip('\n') for line in inClass]
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
outFile.write("site"+'\t'+"codon"+'\t'+"WT_nt"+'\t'+"Variant_nt"+'\t'+"classification"+'\t'+"alert_group"+'\n')


#search for sense sites:


for i in range(0,lenseq-5):
    if seq[i+1:i+6]==sense_rec_site[1:6] and seq[i]!=sense_rec_site[0]:
        for j in range(1,len(inclasslines)):
            inclassthisline = inclasslines[j].split("\t")
            if i == int(inclassthisline[0])-1 and sense_rec_site[0]==inclassthisline[3]:
                outFile.write(inclasslines[j]+'\t'+"group_2"+'\n')

for i in range(0,lenseq-5):
    if seq[i]+seq[i+2:i+6]==sense_rec_site[0]+sense_rec_site[2:6] and seq[i+1]!=sense_rec_site[1]:
        for j in range(1,len(inclasslines)):
            inclassthisline = inclasslines[j].split("\t")
            if i+1 == int(inclassthisline[0])-1 and sense_rec_site[1]==inclassthisline[3]:
                outFile.write(inclasslines[j]+'\t'+"group_2"+'\n')

for i in range(0,lenseq-5):
    if seq[i:i+2]+seq[i+3:i+6]==sense_rec_site[0:2]+sense_rec_site[3:6] and seq[i+2]!=sense_rec_site[2]:
        for j in range(1,len(inclasslines)):
            inclassthisline = inclasslines[j].split("\t")
            if i+2 == int(inclassthisline[0])-1 and sense_rec_site[2]==inclassthisline[3]:
                outFile.write(inclasslines[j]+'\t'+"group_2"+'\n')

for i in range(0,lenseq-5):
    if seq[i:i+3]+seq[i+4:i+6]==sense_rec_site[0:3]+sense_rec_site[4:6] and seq[i+3]!=sense_rec_site[3]:
        for j in range(1,len(inclasslines)):
            inclassthisline = inclasslines[j].split("\t")
            if i+3 == int(inclassthisline[0])-1 and sense_rec_site[3]==inclassthisline[3]:
                outFile.write(inclasslines[j]+'\t'+"group_2"+'\n')

for i in range(0,lenseq-5):
    if seq[i:i+4]+seq[i+5]==sense_rec_site[0:4]+sense_rec_site[5] and seq[i+4]!=sense_rec_site[4]:
        for j in range(1,len(inclasslines)):
            inclassthisline = inclasslines[j].split("\t")
            if i+4 == int(inclassthisline[0])-1 and sense_rec_site[4]==inclassthisline[3]:
                outFile.write(inclasslines[j]+'\t'+"group_2"+'\n')

for i in range(0,lenseq-5):
    if seq[i:i+5]==sense_rec_site[0:5] and seq[i+5]!=sense_rec_site[5]:
        for j in range(1,len(inclasslines)):
            inclassthisline = inclasslines[j].split("\t")
            if i+5 == int(inclassthisline[0])-1 and sense_rec_site[5]==inclassthisline[3]:
                outFile.write(inclasslines[j]+'\t'+"group_2"+'\n')





#search for antisense sites:




for i in range(0,lenseq-5):
    if seq[i+1:i+6]==antisense_rec_site[1:6] and seq[i]!=antisense_rec_site[0]:
        for j in range(1,len(inclasslines)):
            inclassthisline = inclasslines[j].split("\t")
            if i == int(inclassthisline[0])-1 and antisense_rec_site[0]==inclassthisline[3]:
                outFile.write(inclasslines[j]+'\t'+"group_1"+'\n')

for i in range(0,lenseq-5):
    if seq[i]+seq[i+2:i+6]==antisense_rec_site[0]+antisense_rec_site[2:6] and seq[i+1]!=antisense_rec_site[1]:
        for j in range(1,len(inclasslines)):
            inclassthisline = inclasslines[j].split("\t")
            if i+1 == int(inclassthisline[0])-1 and antisense_rec_site[1]==inclassthisline[3]:
                outFile.write(inclasslines[j]+'\t'+"group_1"+'\n')

for i in range(0,lenseq-5):
    if seq[i:i+2]+seq[i+3:i+6]==antisense_rec_site[0:2]+antisense_rec_site[3:6] and seq[i+2]!=antisense_rec_site[2]:
        for j in range(1,len(inclasslines)):
            inclassthisline = inclasslines[j].split("\t")
            if i+2 == int(inclassthisline[0])-1 and antisense_rec_site[2]==inclassthisline[3]:
                outFile.write(inclasslines[j]+'\t'+"group_1"+'\n')

for i in range(0,lenseq-5):
    if seq[i:i+3]+seq[i+4:i+6]==antisense_rec_site[0:3]+antisense_rec_site[4:6] and seq[i+3]!=antisense_rec_site[3]:
        for j in range(1,len(inclasslines)):
            inclassthisline = inclasslines[j].split("\t")
            if i+3 == int(inclassthisline[0])-1 and antisense_rec_site[3]==inclassthisline[3]:
                outFile.write(inclasslines[j]+'\t'+"group_1"+'\n')

for i in range(0,lenseq-5):
    if seq[i:i+4]+seq[i+5]==antisense_rec_site[0:4]+antisense_rec_site[5] and seq[i+4]!=antisense_rec_site[4]:
        for j in range(1,len(inclasslines)):
            inclassthisline = inclasslines[j].split("\t")
            if i+4 == int(inclassthisline[0])-1 and antisense_rec_site[4]==inclassthisline[3]:
                outFile.write(inclasslines[j]+'\t'+"group_1"+'\n')

for i in range(0,lenseq-5):
    if seq[i:i+5]==antisense_rec_site[0:5] and seq[i+5]!=antisense_rec_site[5]:
        for j in range(1,len(inclasslines)):
            inclassthisline = inclasslines[j].split("\t")
            if i+5 == int(inclassthisline[0])-1 and antisense_rec_site[5]==inclassthisline[3]:
                outFile.write(inclasslines[j]+'\t'+"group_1"+'\n')






inFile.close()
inClass.close()
outFile.close()
