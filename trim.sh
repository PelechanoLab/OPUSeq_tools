#!/bin/bash -l
#SBATCH	...

##Import variables from var.txt
maindir=$(grep 'maindir=' var.txt | sed 's/.*maindir\=//')
adap_fw=$(grep 'adap_fw=' var.txt | sed 's/.*adap_fw\=//')
adap_rv=$(grep 'adap_rv=' var.txt | sed 's/.*adap_rv\=//')
adap_err=$(grep 'adap_err=' var.txt | sed 's/.*adap_err\=//')
min_l=$(grep 'min_l=' var.txt | sed 's/.*min_l\=//')

##Load Cutadapt
module load bioinfo-tools cutadapt

cd ${maindir}/ext_fastq

for i in *R1.fastq.gz
do 
SAMPLE=${i%%_R1*}
echo $SAMPLE
date 
R1=${SAMPLE}_R1
R2=${SAMPLE}_R2
R1_trim_fw=${SAMPLE}_trim_fw_R1
R2_trim_fw=${SAMPLE}_trim_fw_R2
R1_trim_rv=${SAMPLE}_trim_R1
R2_trim_rv=${SAMPLE}_trim_R2

cutadapt -e $adap_err -j 8 \
--action=trim \
--minimum-length $min_l \
-g file:$adap_fw \
-G file:$adap_fw \
$R1.fastq.gz \
$R2.fastq.gz \
-o ../trim_fastq/$R1_trim_fw.fastq.gz \
-p ../trim_fastq/$R2_trim_fw.fastq.gz \
&> ../trim_fastq/${SAMPLE}_FW_trim.txt 

cutadapt -e $adap_err -j 8 \
--action=trim \
--minimum-length $min_l \
-a file:$adap_rv \
-A file:$adap_rv \
../trim_fastq/$R1_trim_fw.fastq.gz \
../trim_fastq/$R2_trim_fw.fastq.gz \
-o ../trim_fastq/$R1_trim_rv.fastq.gz \
-p ../trim_fastq/$R2_trim_rv.fastq.gz \
&> ../trim_fastq/${SAMPLE}_RV_trim.txt

rm ../trim_fastq/$R1_trim_fw.fastq.gz 
rm ../trim_fastq/$R2_trim_fw.fastq.gz

done
