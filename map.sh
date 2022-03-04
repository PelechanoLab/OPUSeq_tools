#!/bin/bash -l
#SBATCH	...

##Import variables from var.txt
maindir=$(grep 'maindir=' var.txt | sed 's/.*maindir\=//')
library_type=$(grep 'library_type=' var.txt | sed 's/.*library_type\=//')

##Load BWA, Samtools and iGenomes
module load bioinfo-tools bwa samtools iGenomes

if [[ $library_type == "ds_umi" ]]
then
cd ${maindir}/trim_fastq
	for i in *R1.fastq.gz
	do 
		SAMPLE=${i%%_trim_R1*}
		echo $SAMPLE
		date 
		R1=${SAMPLE}_R1
		R2=${SAMPLE}_R2

		bwa aln -t 8 \
		${IGENOMES_DATA}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa \
		${SAMPLE}_trim_R1.fastq.gz \
		> ../bam/$R1.sai

		bwa aln -t 8 \
		${IGENOMES_DATA}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa \
		${SAMPLE}_trim_R2.fastq.gz \
		> ../bam/$R2.sai

		bwa sampe -a 2000 \
		${IGENOMES_DATA}/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa \
		../bam/$R1.sai \
		../bam/$R2.sai \
		${SAMPLE}_trim_R1.fastq.gz \
		${SAMPLE}_trim_R2.fastq.gz \
		> ../bam/$SAMPLE.sam

		cd ../bam/

		samtools view -@ 8 -b $SAMPLE.sam | samtools sort -@ 8 -o $SAMPLE.bam
		samtools index $SAMPLE.bam $SAMPLE.bam.bai

		rm $SAMPLE.sam $R1.sai $R2.sai

		cd ../trim_fastq/
	done
else
	cd ${maindir}/fastq
	for i in *R1.fastq.gz
	do 
		SAMPLE=${i%%_R1*}
		echo $SAMPLE
		date 
		R1=${SAMPLE}_R1
		R2=${SAMPLE}_R2

		bwa aln -t 8 \
		/sw/data/uppnex/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa \
		$R1.fastq.gz \
		> ../bam/$R1.sai

		bwa aln -t 8 \
		/sw/data/uppnex/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa \
		$R2.fastq.gz \
		> ../bam/$R2.sai

		bwa sampe -a 2000 \
		/sw/data/uppnex/igenomes/Homo_sapiens/UCSC/hg19/Sequence/BWAIndex/genome.fa \
		../bam/$R1.sai \
		../bam/$R2.sai \
		$R1.fastq.gz \
		$R2.fastq.gz \
		> ../bam/$SAMPLE.sam

		cd ../bam/

		samtools view -@ 8 -b $SAMPLE.sam | samtools sort -@ 8 -o $SAMPLE.bam
		samtools index $SAMPLE.bam $SAMPLE.bam.bai

		rm $SAMPLE.sam $R1.sai $R2.sai

		cd ../fastq/
	done
fi

cd ../bam
for i in *.bam
do 
SAMPLE=${i%%.bam*}
ontarget=$(samtools view -@ 8 -c -F 4 -b $i -L ../scr/HRAS_capture_9.bed) 
flanking=$(samtools view -@ 8 -c -F 4 -b $i -L ../scr/HRAS_capture_9_flank.bed)
mapped=$(samtools view -@ 8 -c -F 4 $i)
total=$(samtools view -@ 8 -c $i)

echo total $((total / 2)) >> ${SAMPLE}_stats.txt
echo mapped $((mapped / 2)) >> ${SAMPLE}_stats.txt
echo flanking $((flanking / 2)) >> ${SAMPLE}_stats.txt
echo ontarget $((ontarget / 2)) >> ${SAMPLE}_stats.txt

done





