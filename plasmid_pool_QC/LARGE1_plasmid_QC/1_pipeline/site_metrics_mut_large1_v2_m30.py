# open file

import sys
inFile=open(sys.argv[1],'rt')
outFile=open(sys.argv[2],'wt')
blk=sys.argv[3]

if blk=="1":
    reference="atgctgggaatctgcagggggagacggaaattcttggctgcctcgttgagtcttctctgcatcccagccatcacctggatttacctgttttctgggagcttcgaagatggaaagcccgtgtctctgtcaccgctggagtcccaggcacacagccccaggtacacggcctccagccagcgggagcgcgagagcctggaggtgcgcatgcgcgaggtggaggaggagaa"
if blk=="2":
    reference="ccgcgccctccgcaggcagctcagcctggcccagggccgagccccatcccatcgccgaggcaaccactccaagacctactccatggaggagggcactggagacagcgagaaccttcgggctggcatcgtggcaggcaacagctccgagtgtgggcagcagccggtcgtggagaaatgcgagacaatccacgttgctattgtctgcgccggatacaatgccagccggg"
if blk=="3":
    reference="atgtcgtcaccctggtcaaatccgtcctgttccatagacggaaccctctgcacttccaccttattgctgactccattgcggagcagatcctggccacgctcttccagacctggatggtgcccgctgtgcgtgtggacttctacaatgcagacgagctcaagtctgaagtttcctggatccccaataaacattactctgggatttatggtctgatgaagcttgtcctg"
if blk=="4":
    reference="accaagactcttcctgccaacctggagagagtcatcgtccttgacacggatatcacctttgccactgacattgcagagctgtgggctgtgttccacaagttcaaaggtcagcaagtcctgggcttggtggagaaccagagtgactggtaccttggaaacctgtggaaaaatcaccgcccatggccagcccttggaagaggctacaacacaggggtgatcctgttact"
if blk=="5":
    reference="tctggataagctgcggaagatgaaatgggagcagatgtggaggctgaccgcagagagggagctcatgggcatgctctctacatccttagctgaccaggatattttcaatgccgtcatcaaacaaaaccccttccttgtgtaccagctcccctgcttctggaatgtgcagctgtcagaccacacccgctccgagcagtgctacagagacgtgtctgatctaaaggtca"
if blk=="6":
    reference="ttcactggaactcccccaagaagctccgggtgaagaacaagcatgtggagttttttcgcaacctctacctgaccttcctggagtatgacggcaatcttctgaggcgggaactgtttggctgccccagtgaggctgatgtcaacagtgaaaacctccagaagcagctgtctgagctggacgaggacgacctgtgctatgagttccggcgagagcgcttcactgtccac"
if blk=="7":
    reference="cgcacccacctgtacttcctgcactacgagtatgagcctgcagcagacagcacggacgtcaccctggtcgctcagctgtccatggacaggctccagatgctggaggccatctgcaagcactgggaggggcccatcagcctggccctctacctgtcagacgccgaggcccagcagttcctccgctacgcacagggctctgaggtgcttatgagccgccacaacgtggg"
if blk=="8":
    reference="ctaccacatcgtgtacaaggagggccagttctaccccgtgaacctgctgcgcaacgtggccatgaagcacatcagcactccctacatgttcctgtctgacattgacttcctgcccatgtatgggctctatgagtacctcaggaagtctgtcatccagctcgatcttgccaacaccaagaaagcaatgattgtccccgcgttcgagacactgcgctaccggctgtcct"
if blk=="9":
    reference="tccccaagtcaaaagcggagttgctgtcaatgctggacatggggaccctcttcacattcaggtaccacgtctggacgaaaggccacgcacccacaaacttcgccaagtggcggaccgccaccacgccttaccgggttgagtgggaggccgattttgagccgtatgttgttgtgagacgtgactgcccggagtacgaccggaggtttgtaggctttggctggaacaaa"
if blk=="10":
    reference="gtggctcatatcatggagctggatgtgcaggagtatgagttcattgtgctgcccaacgcctacatgatccacatgcctcatgcccccagcttcgacattaccaagttccgttccaacaagcaataccgcatctgtctcaaaaccctcaaggaagagtttcagcaggacatgtcccgccgctacggctttgctgccctgaaatatctcacagccgagaacaacagc"


#reference sequence



#write a header

outFile.write("chr"+'\t'+"pos"+'\t'+"A"+'\t'+"C"+'\t'+"G"+'\t'+"T"+'\t'+"transition"+'\t'+"transversion"+'\n')

# read the file

SITES= [line.rstrip('\n') for line in inFile]

for i in range(1,len(SITES)):
    NT=SITES[i].split("\t")
    if reference[i-1]=="A" or reference[i-1]=="a":
        NTA=0
        NTC=NT[3]
        NTG=NT[4]
        NTT=NT[5]
        transition=NTG
        transversion=int(NTC)+int(NTT)
    if reference[i-1]=="C" or reference[i-1]=="c":
        NTA=NT[2]
        NTC=0
        NTG=NT[4]
        NTT=NT[5]
        transition=NTT
        transversion=int(NTG)+int(NTA)
    if reference[i-1]=="G" or reference[i-1]=="g":
        NTA=NT[2]
        NTC=NT[3]
        NTG=0
        NTT=NT[5]
        transition=NTA
        transversion=int(NTC)+int(NTT)
    if reference[i-1]=="T" or reference[i-1]=="t":
        NTA=NT[2]
        NTC=NT[3]
        NTG=NT[4]
        NTT=0
        transition=NTC
        transversion=int(NTA)+int(NTG)

    outFile.write(NT[0]+'\t'+str(int(NT[1])-30)+'\t'+str(NTA)+'\t'+str(NTC)+'\t'+str(NTG)+'\t'+str(NTT)+'\t'+str(transition)+'\t'+str(transversion)+'\n')

inFile.close()
outFile.close()
