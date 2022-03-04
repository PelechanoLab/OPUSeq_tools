#!/bin/bash -l
#SBATCH	...

repeat(){
	local start=1
	local end=${1:-80}
	local str="${2:-=}"
	local range=$(seq $start $end)
	for i in $range ; do echo -n "${str}"; done
}

##Import variables from var.txt
maindir=$(grep 'maindir=' var.txt | sed 's/.*maindir\=//')
umi_l=$(grep 'umi_l=' var.txt | sed 's/.*umi_l\=//')
umi_pattern=$(repeat $umi_l 'N'; echo) ##the pattern for UMI-tools extract to look for: N repeated $umi_l times.

##Load UMI tools 
module load bioinfo-tools umi_tools 

cd ${maindir}/fastq

for i in *R1.fastq.gz
do 
SAMPLE=${i%%_R1*}
R1_in=${SAMPLE}_R1.fastq.gz
R2_in=${SAMPLE}_R2.fastq.gz
R1_out=${SAMPLE}_R1.fastq.gz
R2_out=${SAMPLE}_R2.fastq.gz
echo $SAMPLE
date
umi_tools extract \
--extract-method=string \
--bc-pattern=$umi_pattern \
--bc-pattern2=$umi_pattern \
-I $R1_in \
-S ../ext_fastq/$R1_out \
--read2-in=$R2_in \
--read2-out=../ext_fastq/$R2_out \
-L ../ext_fastq/$SAMPLE.log 
done






