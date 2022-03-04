#!/bin/bash -l
#SBATCH	..

##Import variables from var.txt
maindir=$(grep 'maindir=' var.txt | sed 's/.*maindir\=//')
umi_cons=$(grep 'umi_cons=' var.txt | sed 's/.*umi_cons\=//')
bound_only_cons=$(grep 'bound_only_cons=' var.txt | sed 's/.*bound_only_cons\=//')
minmem=$(grep 'minmem=' var.txt | sed 's/.*minmem\=//')
dcs=$(grep 'dcs=' var.txt | sed 's/.*dcs\=//')
sscs=$(grep 'sscs=' var.txt | sed 's/.*sscs\=//')

cd $maindir/pp_bam
list_bam=$(ls *.bam)

##Load Python 3 and Biopython.
module load bioinfo-tools python biopython

if $umi_cons
then
	for i in $list_bam
	do
		SAMPLE=${i%%_PP.bam*}
		echo $i 
		date
		if $dcs
		then
			if $sscs
			then
				python ../scr/consensus_read_mod.py --input $i --UMI_in_name --prefix $SAMPLE --minmem $minmem --numCores 2 --tagstats --write-sscs
			else
				python ../scr/consensus_read_mod.py --input $i --UMI_in_name --prefix $SAMPLE --minmem $minmem --numCores 2 --tagstats
			fi
		elif $sscs
		then
			python ../scr/consensus_read_mod.py --input $i --UMI_in_name --prefix $SAMPLE --minmem $minmem --numCores 2 --tagstats --without-dcs --write-sscs
		fi
	done
fi 

if $bound_only_cons
then
	for i in $list_bam
	do
		SAMPLE=${i%%_PP.bam*}
		echo $i 
		date
		if $dcs
		then
			if $sscs
			then
				python ../scr/consensus_read_mod.py --input $i --prefix ${SAMPLE}_noUMI --minmem $minmem --numCores 2 --tagstats --write-sscs
			else
				python ../scr/consensus_read_mod.py --input $i --prefix ${SAMPLE}_noUMI --minmem $minmem --numCores 2 --tagstats
			fi
		elif $sscs
		then
			python ../scr/consensus_read_mod.py --input $i --prefix ${SAMPLE}_noUMI --minmem $minmem --numCores 2 --tagstats --without-dcs --write-sscs
		fi
	done
fi

##Move the resulting fastq and stats files to the appropriate folders:

if $dcs
then

	for i in *dcs.fq.gz
	do
	mv $i ../dcs_fastq/$i
	done

fi 

if $sscs
then

	for i in *sscs.fq.gz
	do
	mv $i ../sscs_fastq/$i
	done

fi

for i in *.txt
do 
	mv $i ../cmstats/$i
done

for i in *.png
do 
	mv $i ../cmstats/$i
done
