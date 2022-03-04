#!/bin/bash -l
#SBATCH ...

##Import variables from var.txt
maindir=$(grep 'maindir=' var.txt | sed 's/.*maindir\=//')
edit_dist=$(grep 'edit_dist=' var.txt | sed 's/.*edit_dist\=//')

##Load UMI-tools 
module load bioinfo-tools umi_tools

cd $maindir/bam

for i in *.bam
do
	SAMPLE=${i%%.bam*}
	echo $i 
	date
	umi_tools group -I $i \
	--paired --edit-distance-threshold $edit_dist \
	--chimeric-pairs=discard \
	--unpaired-reads=discard \
	--group-out=../tsv/${SAMPLE}.tsv \
	--log=../tsv/${SAMPLE}.log 
done


