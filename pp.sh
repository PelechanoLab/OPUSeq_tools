#!/bin/bash -l
#SBATCH	...

##Import variables from var.txt
maindir=$(grep 'maindir=' var.txt | sed 's/.*maindir\=//')
library_type=$(grep 'library_type=' var.txt | sed 's/.*library_type\=//')
pp_mapq=$(grep 'pp_mapq=' var.txt | sed 's/.*pp_mapq\=//')

cd $maindir/bam
list_bam=$(ls *.bam)

##Load python 3 and Biopython
module load bioinfo-tools python biopython

if [[ $library_type == "ds_umi" ]]
then
	for i in $list_bam
		do
			SAMPLE=${i%%.bam*}
			echo $i 
			date
			python ../scr/correct_pair.py --input $i --out ../pp_bam/${SAMPLE}_PP.bam --stat ../pp_bam/${SAMPLE}_stats.txt --Q $pp_mapq --UMI-group ../tsv/${SAMPLE}.tsv 
		done
else
	for i in $list_bam
		do
			SAMPLE=${i%%.bam*}
			echo $i 
			date
			python ../scr/correct_pair.py --input $i --out ../pp_bam/${SAMPLE}_PP.bam --stat ../pp_bam/${SAMPLE}_stats.txt --Q $pp_mapq 
		done
fi


