#!/bin/bash -l
#SBATCH	...

##Import variables from var.txt
maindir=$(grep 'maindir=' var.txt | sed 's/.*maindir\=//')
dcs=$(grep 'dcs=' var.txt | sed 's/.*dcs\=//')
sscs=$(grep 'sscs=' var.txt | sed 's/.*sscs\=//')

##Load Bowtie2, Samtools, and iGenomes.
module load bioinfo-tools bowtie2 samtools iGenomes

if $sscs
then
	cd $maindir/sscs_fastq
	for i in *read1*fq.gz
	do
	SAMPLE=${i%%_read*_sscs.fq.gz*}
	echo $SAMPLE
	date

	R1=${SAMPLE}_read1_sscs.fq.gz
	R2=${SAMPLE}_read2_sscs.fq.gz

	bowtie2 -p 4 --local -x ${IGENOMES_DATA}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome \
	--np 0 --n-ceil L,0,1 -1 $R1 -2 $R2 -S ../sscs_bam/${SAMPLE}_SSCS.sam 

	cd ../sscs_bam/
	samtools view -@ 4 -b ${SAMPLE}_SSCS.sam | samtools sort -@ 4 -o ${SAMPLE}_SSCS.bam
	samtools index ${SAMPLE}_SSCS.bam ${SAMPLE}_SSCS.bam.bai
	rm ${SAMPLE}_SSCS.sam 

	cd ../sscs_fastq
	done
fi


if $dcs
then
	cd $maindir/dcs_fastq
	for i in *read1*fq.gz
	do
	SAMPLE=${i%%_read*_dcs.fq.gz*}
	echo $SAMPLE
	date

	R1=${SAMPLE}_read1_dcs.fq.gz
	R2=${SAMPLE}_read2_dcs.fq.gz

	bowtie2 -p 4 --local -x ${IGENOMES_DATA}/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome \
	--np 0 --n-ceil L,0,1 -1 $R1 -2 $R2 -S ../dcs_bam/${SAMPLE}_DCS.sam 

	cd ../dcs_bam/
	samtools view -@ 4 -b ${SAMPLE}_DCS.sam | samtools sort -@ 4 -o ${SAMPLE}_DCS.bam
	samtools index ${SAMPLE}_DCS.bam ${SAMPLE}_DCS.bam.bai
	rm ${SAMPLE}_DCS.sam 

	cd ../dcs_fastq
	done
fi
