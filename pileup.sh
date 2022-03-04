#!/bin/bash -l
#SBATCH	...

##Import variables from var.txt
maindir=$(grep 'maindir=' var.txt | sed 's/.*maindir\=//')
dcs=$(grep 'dcs=' var.txt | sed 's/.*dcs\=//')
sscs=$(grep 'sscs=' var.txt | sed 's/.*sscs\=//')
raw=$(grep 'raw=' var.txt | sed 's/.*raw\=//')
pileup_q=$(grep 'pileup_q=' var.txt | sed 's/.*pileup_q\=//')
bam_trim=$(grep 'bam_trim=' var.txt | sed 's/.*bam_trim\=//')

##Load Samtools and Python 3
module load bioinfo-tools samtools python

##Perform pileup with samtools. Do not filter DCS and SSCS bam files by Phred quality, only the raw files.

if $dcs
then
	cd $maindir/dcs_bam
	for i in *.bam
	do
		SAMPLE=${i%%.bam*}
		samtools mpileup -f /sw/data/uppnex/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
		--max-depth 0 --count-orphans -q 1 -B $i > ../pileup/${SAMPLE}_pileup.txt
	done
fi

if $sscs
then
	cd $maindir/sscs_bam
	for i in *.bam
	do
		SAMPLE=${i%%.bam*}
		samtools mpileup -f /sw/data/uppnex/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
		--max-depth 0 --count-orphans -q 1 -B $i > ../pileup/${SAMPLE}_pileup.txt
	done
fi

if $raw
then
	cd $maindir/pp_bam
	for i in *.bam
	do
		SAMPLE=${i%%.bam*}
		samtools mpileup -f /sw/data/uppnex/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
		--max-depth 0 --count-orphans -q 1 -Q $pileup_q -B $i > ../pileup/${SAMPLE}_raw_pileup.txt
	done
fi

##Collapse the pileups into .csv tables; delete original pileups to save space
cd $maindir/pileup
for i in *pileup.txt
do
	SAMPLE=${i%%_pileup.txt*}
	python ../scr/collapse_pileup.py --input $i --output ${SAMPLE}_vars.csv >> ${SAMPLE}.log
	rm $i
done

##If doing bam trimming, perform the same steps for trimmed bam files:

if $bam_trim
then
mkdir $maindir/pileup_trim
if $dcs
then
cd $maindir/dcs_bam_trim
for i in *.bam
do
	SAMPLE=${i%%.bam*}
	samtools mpileup -f /sw/data/uppnex/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
	--max-depth 0 --count-orphans -q 1 -B $i > ../pileup_trim/${SAMPLE}_pileup.txt
done
fi

if $sscs
then
cd $maindir/sscs_bam_trim
for i in *.bam
do
SAMPLE=${i%%.bam*}
samtools mpileup -f /sw/data/uppnex/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
--max-depth 0 --count-orphans -q 1 -B $i > ../pileup_trim/${SAMPLE}_pileup.txt
done
fi

if $raw
then
cd $maindir/pp_bam_trim
for i in *.bam
do
	SAMPLE=${i%%.bam*}
	samtools mpileup -f /sw/data/uppnex/igenomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa \
	--max-depth 0 --count-orphans -q 1 -Q $pileup_q -B $i > ../pileup_trim/${SAMPLE}_raw_pileup.txt
done
fi

cd $maindir/pileup_trim
for i in *pileup.txt
do
SAMPLE=${i%%_pileup.txt*}
python ../scr/collapse_pileup.py --input $i --output ${SAMPLE}_vars.csv >> ${SAMPLE}.log
rm $i
done

fi







