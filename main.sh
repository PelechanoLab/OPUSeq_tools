#!/bin/bash -l
#SBATCH	...

##All the variables must be specified manually in a txt file called "var.txt"
##Get the variables needed for this script:

maindir=$(grep 'maindir' var.txt | sed 's/.*maindir\=//')
library_type=$(grep 'library_type' var.txt | sed 's/.*library_type\=//')
dcs=$(grep 'dcs' var.txt | sed 's/.*dcs\=//')
sscs=$(grep 'sscs' var.txt | sed 's/.*sscs\=//')
umi_cons=$(grep 'umi_cons' var.txt | sed 's/.*umi_cons\=//')
bound_only_cons=$(grep 'bound_only_cons' var.txt | sed 's/.*bound_only_cons\=//')
raw=$(grep 'raw' var.txt | sed 's/.*raw\=//')

##Make the necessary directories.

cd $maindir

##First the general directories, needed for all types of analysis: (a fastq folder containing the fastq files must already be created in main)
mkdir bam pp_bam pileup

##Then the folders general to OPUSeq:
if [[ $library_type == "ds_umi" ]]
then
	mkdir ext_fastq trim_fastq
fi

##If making consensus, make cmstats (consensus making stats) folder:
if $dcs || $sscs
then 
	mkdir cmstats
fi

##Make tsv dir for all cases where UMI is used:
if $umi_cons
then
	mkdir tsv
fi 

##Make dcs- and sscs-specific folders for analysis with and/or without UMIs:

if $dcs
then
	mkdir dcs_fastq dcs_bam
fi

if $sscs
then 
	mkdir sscs_fastq sscs_bam
fi

##If library_type = no_umi, then no analysis using UMIs can be performed.
##Therefore, if library_type = no_umi and umi_cons = T, raise error
if [[ $library_type == "no_umi" ]] && $umi_cons 
then
	echo 'If the library contains no UMIs, then analysis using UMIs cannot be performed.
	I.e., $library_type cannot be no_umi at the same time as $umi_cons is set to TRUE.'
	exit 1
fi

##If at least one of dcs and sscs is set to true, then at least one of 
##umi_cons or bound_only_cons must be set to true. 
if $dcs || $sscs
then
if ! $umi_cons && ! $bound_only_cons
then
echo 'If at least one of $dcs and $sscs is set to TRUE, then at least one of 
$umi_cons or $bound_only_cons must be set to TRUE.'
exit 1
fi
fi

#####################################

cd scr

if [[ $library_type == "ds_umi" ]]
then
	ID1=$(sbatch umi_extract.sh | grep -P '[0-9]*' -o)
	ID2=$(sbatch --dependency=afterok:$ID1 trim.sh | grep -P '[0-9]*' -o)
	ID3=$(sbatch --dependency=afterok:$ID2 map.sh | grep -P '[0-9]*' -o)
	ID4=$(sbatch --dependency=afterok:$ID3 group.sh | grep -P '[0-9]*' -o)
	ID5=$(sbatch --dependency=afterok:$ID4 pp.sh | grep -P '[0-9]*' -o)
	if $dcs || $sscs
	then
		ID6=$(sbatch --dependency=afterok:$ID5 consensus.sh | grep -P '[0-9]*' -o)
		ID7=$(sbatch --dependency=afterok:$ID6 consensus_map.sh | grep -P '[0-9]*' -o)
		ID8=$(sbatch --dependency=afterok:$ID7 bam_filter_trim.sh | grep -P '[0-9]*' -o) 
		sbatch --dependency=afterok:$ID8 pileup.sh
	elif $raw
	then
		ID8=$(sbatch --dependency=afterok:$ID5 bam_filter_trim.sh | grep -P '[0-9]*' -o) 
		sbatch --dependency=afterok:$ID8 pileup.sh
	fi
else
	ID1=$(sbatch map.sh | grep -P '[0-9]*' -o)
	ID2=$(sbatch --dependency=afterok:$ID1 pp.sh | grep -P '[0-9]*' -o)
	if $dcs || $sscs
	then 
		ID3=$(sbatch --dependency=afterok:$ID2 consensus.sh | grep -P '[0-9]*' -o)
		ID4=$(sbatch --dependency=afterok:$ID3 consensus_map.sh | grep -P '[0-9]*' -o)
		ID5=$(sbatch --dependency=afterok:$ID4 bam_filter_trim.sh | grep -P '[0-9]*' -o)
		sbatch --dependency=afterok:$ID5 pileup.sh
	elif $raw
	then
		ID5=$(sbatch --dependency=afterok:$ID2 bam_filter_trim.sh | grep -P '[0-9]*' -o)
		sbatch --dependency=afterok:$ID5 pileup.sh
	fi
fi





