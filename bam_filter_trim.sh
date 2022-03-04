#!/bin/bash -l
#SBATCH	...

maindir=$(grep 'maindir=' var.txt | sed 's/.*maindir\=//')
dcs=$(grep 'dcs=' var.txt | sed 's/.*dcs\=//')
sscs=$(grep 'sscs=' var.txt | sed 's/.*sscs\=//')
bam_trim=$(grep 'bam_trim=' var.txt | sed 's/.*bam_trim\=//')
bam_trim_bp=$(grep 'bam_trim_bp=' var.txt | sed 's/.*bam_trim_bp\=//')
bamutil_path=$(grep 'bamutil_path=' var.txt | sed 's/.*bamutil_path\=//')

##Load Samtools
module load bioinfo-tools samtools 

if $dcs 
then
##Filter the bam files
cd $maindir/dcs_bam
list_bam=$(ls *.bam)
for i in $list_bam
do
SAMPLE=${i%%.bam*}
samtools view -@ 8 -q 36 -b $i -o ${SAMPLE}_filter.bam
rm $i
mv ${SAMPLE}_filter.bam $i
samtools index $i $i.bai
done
if $bam_trim
then
##Trim bases from each end of each read in the bam files using bamUtil
mkdir $maindir/dcs_bam_trim
cd $maindir/dcs_bam
list_bam=$(ls *.bam)
cd $bamutil_path
for i in $list_bam
do
./bam trimBam $maindir/dcs_bam/$i $maindir/dcs_bam_trim/$i $bam_trim_bp
done
##Index the resulting trimmed files
cd $maindir/dcs_bam_trim
list_bam=$(ls *.bam)
for i in $list_bam
do
samtools index -@ 8 $i $i.bai
done
fi
fi

##Repeat for sscs and raw if those variables are set to TRUE

if $sscs
then
cd $maindir/sscs_bam
list_bam=$(ls *.bam)
for i in $list_bam
do
SAMPLE=${i%%.bam*}
samtools view -@ 8 -q 36 -b $i -o ${SAMPLE}_filter.bam
rm $i
mv ${SAMPLE}_filter.bam $i
samtools index $i $i.bai
done

if $bam_trim
then
mkdir $maindir/sscs_bam_trim
cd $maindir/sscs_bam
list_bam=$(ls *.bam)
cd $bamutil_path
for i in $list_bam
do
./bam trimBam $maindir/sscs_bam/$i $maindir/sscs_bam_trim/$i $bam_trim_bp
done

cd $maindir/sscs_bam_trim
list_bam=$(ls *.bam)
for i in $list_bam
do
samtools index -@ 8 $i $i.bai
done
fi
fi

##Do not filter raw files, since they are already filtered in pp.sh

if $raw
then

if $bam_trim
mkdir $maindir/pp_bam_trim
then
cd $maindir/pp_bam
list_bam=$(ls *.bam)
cd $bamutil_path
for i in $list_bam
do
./bam trimBam $maindir/pp_bam/$i $maindir/pp_bam_trim/$i $bam_trim_bp
done

cd $maindir/pp_bam_trim
list_bam=$(ls *.bam)
for i in $list_bam
do
samtools index -@ 8 $i $i.bai
done
fi
fi














