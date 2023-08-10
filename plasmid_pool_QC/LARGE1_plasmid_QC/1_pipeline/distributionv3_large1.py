# open file

import sys
inFile=open(sys.argv[1],'rt')
inFiletotal=open(sys.argv[2],'rt')
outFile=open(sys.argv[3],'wt')
blk=sys.argv[4]

# set the counter

readbin0=-227

if blk=="10":
    readbin0=-225
    

readbin1=0
readbin2=0
readbin3=0
readbin4=0
readbin5=0




# read the file
NUMBERStotal= [line.rstrip('\n') for line in inFiletotal]
NUMBERS= [line.rstrip('\n') for line in inFile]
for i in range(1,len(NUMBERS)):
    NUMBER=NUMBERS[i].split("\t")
    NUMBERtotal=NUMBERStotal[i].split("\t")
    ttreads=int(NUMBERtotal[2])+int(NUMBERtotal[3])+int(NUMBERtotal[4])+int(NUMBERtotal[5])

    for j in range(2,6):
        if float(NUMBER[j])/float(ttreads)==0:
            readbin0+=1

        if float(NUMBER[j])/float(ttreads)>0 and float(NUMBER[j])/float(ttreads)<=0.0005:
            readbin1+=1

        if float(NUMBER[j])/float(ttreads)>0.0005 and float(NUMBER[j])/float(ttreads)<=0.001:
            readbin2+=1

        if float(NUMBER[j])/float(ttreads)>0.001 and float(NUMBER[j])/float(ttreads)<=0.0015:
            readbin3+=1

        if float(NUMBER[j])/float(ttreads)>0.0015 and float(NUMBER[j])/float(ttreads)<=0.002:
            readbin4+=1

        if float(NUMBER[j])/float(ttreads)>0.002:
            readbin5+=1

outFile.write("readbin0"+'\t'+str(readbin0)+'\n'+"readbin1"+'\t'+str(readbin1)+'\n'+"readbin2"+'\t'+str(readbin2)+'\n'+"readbin3"+'\t'+str(readbin3)+'\n'+"readbin4"+'\t'+str(readbin4)+'\n'+"readbin5"+'\t'+str(readbin5)+'\n')

inFile.close()
inFiletotal.close()
outFile.close()
