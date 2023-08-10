#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=QCLARGE1
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=64000
#SBATCH --time=36:00:00


# change to the working directory:
cd /gpfs/ycga/scratch60/lek/km2286/QCLARGE1;

# source my bashrc file:
source ~/.bashrc

# run the commands



# build a CDS reference
module load BWA/0.7.15-foss-2016a

bwa index LARGE1_CDS_allcap_30s30e.fasta

module load SAMtools

samtools faidx LARGE1_CDS_allcap_30s30e.fasta





# demultiplexing
module load FASTX-Toolkit/0.0.14-GCCcore-10.2.0
cat large1plaspoo_R1_001.fastq | /ycga-gpfs/apps/hpc/software/FASTX-Toolkit/0.0.14-GCCcore-10.2.0/bin/fastx_barcode_splitter.pl --bcfile barcode --bol --mismatches 5 --prefix R1.  --suffix .fastq

# edit file name in find_R2

chmod 777 find_R2.pl
chmod 777 make_all_R2.pl
./make_all_R2.pl




# alignment


module load BWA/0.7.15-foss-2016a

bwa mem -M -t 8 -R "@RG\tID:large1\tPL:ILLUMINA\tSM:large1" LARGE1_CDS_allcap_30s30e.fasta R1.no1.fastq R2.no1.fastq > no1.sam 2> no1.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:large1\tPL:ILLUMINA\tSM:large1" LARGE1_CDS_allcap_30s30e.fasta R1.no2.fastq R2.no2.fastq > no2.sam 2> no2.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:large1\tPL:ILLUMINA\tSM:large1" LARGE1_CDS_allcap_30s30e.fasta R1.no3.fastq R2.no3.fastq > no3.sam 2> no3.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:large1\tPL:ILLUMINA\tSM:large1" LARGE1_CDS_allcap_30s30e.fasta R1.no4.fastq R2.no4.fastq > no4.sam 2> no4.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:large1\tPL:ILLUMINA\tSM:large1" LARGE1_CDS_allcap_30s30e.fasta R1.no5.fastq R2.no5.fastq > no5.sam 2> no5.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:large1\tPL:ILLUMINA\tSM:large1" LARGE1_CDS_allcap_30s30e.fasta R1.no6.fastq R2.no6.fastq > no6.sam 2> no6.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:large1\tPL:ILLUMINA\tSM:large1" LARGE1_CDS_allcap_30s30e.fasta R1.no7.fastq R2.no7.fastq > no7.sam 2> no7.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:large1\tPL:ILLUMINA\tSM:large1" LARGE1_CDS_allcap_30s30e.fasta R1.no8.fastq R2.no8.fastq > no8.sam 2> no8.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:large1\tPL:ILLUMINA\tSM:large1" LARGE1_CDS_allcap_30s30e.fasta R1.no9.fastq R2.no9.fastq > no9.sam 2> no9.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:large1\tPL:ILLUMINA\tSM:large1" LARGE1_CDS_allcap_30s30e.fasta R1.no10.fastq R2.no10.fastq > no10.sam 2> no10.bwa.sterr.log

# make bam files

module load picard/2.9.0-Java-1.8.0_121

java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=no1.sam OUTPUT= no1.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=no2.sam OUTPUT= no2.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=no3.sam OUTPUT= no3.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=no4.sam OUTPUT= no4.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=no5.sam OUTPUT= no5.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=no6.sam OUTPUT= no6.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=no7.sam OUTPUT= no7.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=no8.sam OUTPUT= no8.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=no9.sam OUTPUT= no9.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=no10.sam OUTPUT= no10.bam SORT_ORDER=coordinate

module load SAMtools

samtools index no1.bam
samtools index no2.bam
samtools index no3.bam
samtools index no4.bam
samtools index no5.bam
samtools index no6.bam
samtools index no7.bam
samtools index no8.bam
samtools index no9.bam
samtools index no10.bam

# do the site_metrics

chmod 777 site_metrics

./site_metrics -r LARGE1_CDS:31-257 -s -d 0 -q 37 no1.bam> no1.out1
./site_metrics -r LARGE1_CDS:258-484 -s -d 0 -q 37 no2.bam> no2.out1
./site_metrics -r LARGE1_CDS:485-711 -s -d 0 -q 37 no3.bam> no3.out1
./site_metrics -r LARGE1_CDS:712-938 -s -d 0 -q 37 no4.bam> no4.out1
./site_metrics -r LARGE1_CDS:939-1165 -s -d 0 -q 37 no5.bam> no5.out1
./site_metrics -r LARGE1_CDS:1166-1392 -s -d 0 -q 37 no6.bam> no6.out1
./site_metrics -r LARGE1_CDS:1393-1619 -s -d 0 -q 37 no7.bam> no7.out1
./site_metrics -r LARGE1_CDS:1620-1846 -s -d 0 -q 37 no8.bam> no8.out1
./site_metrics -r LARGE1_CDS:1847-2073 -s -d 0 -q 37 no9.bam> no9.out1
./site_metrics -r LARGE1_CDS:2074-2298 -s -d 0 -q 37 no10.bam> no10.out1



# distribution analysis

python site_metrics_mut_large1_v2_m30.py no1.out1 no1.out2 1
python site_metrics_mut_large1_v2_m30.py no2.out1 no2.out2 2
python site_metrics_mut_large1_v2_m30.py no3.out1 no3.out2 3
python site_metrics_mut_large1_v2_m30.py no4.out1 no4.out2 4
python site_metrics_mut_large1_v2_m30.py no5.out1 no5.out2 5
python site_metrics_mut_large1_v2_m30.py no6.out1 no6.out2 6
python site_metrics_mut_large1_v2_m30.py no7.out1 no7.out2 7
python site_metrics_mut_large1_v2_m30.py no8.out1 no8.out2 8
python site_metrics_mut_large1_v2_m30.py no9.out1 no9.out2 9
python site_metrics_mut_large1_v2_m30.py no10.out1 no10.out2 10

python distributionv3_large1.py no1.out2 no1.out1 no1out 1
python distributionv3_large1.py no2.out2 no2.out1 no2out 2
python distributionv3_large1.py no3.out2 no3.out1 no3out 3
python distributionv3_large1.py no4.out2 no4.out1 no4out 4
python distributionv3_large1.py no5.out2 no5.out1 no5out 5
python distributionv3_large1.py no6.out2 no6.out1 no6out 6
python distributionv3_large1.py no7.out2 no7.out1 no7out 7
python distributionv3_large1.py no8.out2 no8.out1 no8out 8
python distributionv3_large1.py no9.out2 no9.out1 no9out 9
python distributionv3_large1.py no10.out2 no10.out1 no10out 10
