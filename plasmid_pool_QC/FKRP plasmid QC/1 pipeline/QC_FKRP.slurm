#!/bin/bash
#SBATCH --partition=general
#SBATCH --job-name=QCFKRP
#SBATCH --ntasks=1 --nodes=1
#SBATCH --mem-per-cpu=64000
#SBATCH --time=36:00:00


# change to the working directory:
cd /gpfs/ycga/scratch60/lek/km2286/QCFKRP;

# source my bashrc file:
source ~/.bashrc

# run the commands



# build a CDS reference
module load BWA/0.7.15-foss-2016a

bwa index FKRPforREF.fasta

module load SAMtools

samtools faidx FKRPforREF.fasta





# demultiplexing
module load FASTX-Toolkit/0.0.14-GCCcore-10.2.0

cat 2pool3pool_R1_001.fastq| /ycga-gpfs/apps/hpc/software/FASTX-Toolkit/0.0.14-GCCcore-10.2.0/bin/fastx_barcode_splitter.pl --bcfile barcode23 --bol --mismatches 5 --prefix R1.  --suffix .fastq
cat sample1_R1_001.fastq| /ycga-gpfs/apps/hpc/software/FASTX-Toolkit/0.0.14-GCCcore-10.2.0/bin/fastx_barcode_splitter.pl --bcfile barcode14 --bol --mismatches 5 --prefix R1.  --suffix .fastq
cat sample2_R1_001.fastq| /ycga-gpfs/apps/hpc/software/FASTX-Toolkit/0.0.14-GCCcore-10.2.0/bin/fastx_barcode_splitter.pl --bcfile barcode56 --bol --mismatches 5 --prefix R1.  --suffix .fastq







# edit file name in find_R2

chmod 777 find_R2_14.pl
chmod 777 find_R2_23.pl
chmod 777 find_R2_56.pl

chmod 777 make_all_R2_14.pl
chmod 777 make_all_R2_23.pl
chmod 777 make_all_R2_56.pl

./make_all_R2_14.pl
./make_all_R2_23.pl
./make_all_R2_56.pl




# alignment


module load BWA/0.7.15-foss-2016a

bwa mem -M -t 8 -R "@RG\tID:FKRP\tPL:ILLUMINA\tSM:FKRP" FKRPforREF.fasta R1.14_Blk1f.fastq R2.14_Blk1f.fastq > blk1f.sam 2> blk1f.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:FKRP\tPL:ILLUMINA\tSM:FKRP" FKRPforREF.fasta R1.14_Blk1r.fastq R2.14_Blk1r.fastq > blk1r.sam 2> blk1r.bwa.sterr.log

bwa mem -M -t 8 -R "@RG\tID:FKRP\tPL:ILLUMINA\tSM:FKRP" FKRPforREF.fasta R1.23_Blk2f.fastq R2.23_Blk2f.fastq > blk2f.sam 2> blk2f.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:FKRP\tPL:ILLUMINA\tSM:FKRP" FKRPforREF.fasta R1.23_Blk2r.fastq R2.23_Blk2r.fastq > blk2r.sam 2> blk2r.bwa.sterr.log

bwa mem -M -t 8 -R "@RG\tID:FKRP\tPL:ILLUMINA\tSM:FKRP" FKRPforREF.fasta R1.23_Blk3f.fastq R2.23_Blk3f.fastq > blk3f.sam 2> blk3f.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:FKRP\tPL:ILLUMINA\tSM:FKRP" FKRPforREF.fasta R1.23_Blk3r.fastq R2.23_Blk3r.fastq > blk3r.sam 2> blk3r.bwa.sterr.log

bwa mem -M -t 8 -R "@RG\tID:FKRP\tPL:ILLUMINA\tSM:FKRP" FKRPforREF.fasta R1.14_Blk4f.fastq R2.14_Blk4f.fastq > blk4f.sam 2> blk4f.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:FKRP\tPL:ILLUMINA\tSM:FKRP" FKRPforREF.fasta R1.14_Blk4r.fastq R2.14_Blk4r.fastq > blk4r.sam 2> blk4r.bwa.sterr.log

bwa mem -M -t 8 -R "@RG\tID:FKRP\tPL:ILLUMINA\tSM:FKRP" FKRPforREF.fasta R1.56_Blk5f.fastq R2.56_Blk5f.fastq > blk5f.sam 2> blk5f.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:FKRP\tPL:ILLUMINA\tSM:FKRP" FKRPforREF.fasta R1.56_Blk5r.fastq R2.56_Blk5r.fastq > blk5r.sam 2> blk5r.bwa.sterr.log

bwa mem -M -t 8 -R "@RG\tID:FKRP\tPL:ILLUMINA\tSM:FKRP" FKRPforREF.fasta R1.56_Blk6f.fastq R2.56_Blk6f.fastq > blk6f.sam 2> blk6f.bwa.sterr.log
bwa mem -M -t 8 -R "@RG\tID:FKRP\tPL:ILLUMINA\tSM:FKRP" FKRPforREF.fasta R1.56_Blk6r.fastq R2.56_Blk6r.fastq > blk6r.sam 2> blk6r.bwa.sterr.log


sed -i -e 1,3d blk1r.sam
sed -i -e 1,3d blk2r.sam
sed -i -e 1,3d blk3r.sam
sed -i -e 1,3d blk4r.sam
sed -i -e 1,3d blk5r.sam
sed -i -e 1,3d blk6r.sam

cat blk1f.sam blk1r.sam > blk1.sam
cat blk2f.sam blk2r.sam > blk2.sam
cat blk3f.sam blk3r.sam > blk3.sam
cat blk4f.sam blk4r.sam > blk4.sam
cat blk5f.sam blk5r.sam > blk5.sam
cat blk6f.sam blk6r.sam > blk6.sam


# make bam files

module load picard/2.9.0-Java-1.8.0_121

java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=blk1.sam OUTPUT= blk1.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=blk2.sam OUTPUT= blk2.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=blk3.sam OUTPUT= blk3.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=blk4.sam OUTPUT= blk4.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=blk5.sam OUTPUT= blk5.bam SORT_ORDER=coordinate
java -Xmx14G -jar $EBROOTPICARD/picard.jar SortSam INPUT=blk6.sam OUTPUT= blk6.bam SORT_ORDER=coordinate


module load SAMtools

samtools index blk1.bam
samtools index blk2.bam
samtools index blk3.bam
samtools index blk4.bam
samtools index blk5.bam
samtools index blk6.bam



# do the site_metrics

chmod 777 site_metrics

./site_metrics -r FKRP_1to1485_s30e30:31-280 -s -d 0 -q 37 blk1.bam> blk1.out1
./site_metrics -r FKRP_1to1485_s30e30:281-530 -s -d 0 -q 37 blk2.bam> blk2.out1
./site_metrics -r FKRP_1to1485_s30e30:531-780 -s -d 0 -q 37 blk3.bam> blk3.out1
./site_metrics -r FKRP_1to1485_s30e30:781-1030 -s -d 0 -q 37 blk4.bam> blk4.out1
./site_metrics -r FKRP_1to1485_s30e30:1031-1280 -s -d 0 -q 37 blk5.bam> blk5.out1
./site_metrics -r FKRP_1to1485_s30e30:1281-1515 -s -d 0 -q 37 blk6.bam> blk6.out1



# distribution analysis

python site_metrics_mut_FKRP_v2_m30.py blk1.out1 blk1.out2 1
python site_metrics_mut_FKRP_v2_m30.py blk2.out1 blk2.out2 2
python site_metrics_mut_FKRP_v2_m30.py blk3.out1 blk3.out2 3
python site_metrics_mut_FKRP_v2_m30.py blk4.out1 blk4.out2 4
python site_metrics_mut_FKRP_v2_m30.py blk5.out1 blk5.out2 5
python site_metrics_mut_FKRP_v2_m30.py blk6.out1 blk6.out2 6




python distributionv3_FKRP.py blk1.out2 blk1.out1 blk1out 1
python distributionv3_FKRP.py blk2.out2 blk2.out1 blk2out 2
python distributionv3_FKRP.py blk3.out2 blk3.out1 blk3out 3
python distributionv3_FKRP.py blk4.out2 blk4.out1 blk4out 4
python distributionv3_FKRP.py blk5.out2 blk5.out1 blk5out 5
python distributionv3_FKRP.py blk6.out2 blk6.out1 blk6out 6
