# OPUSeq_tools

## Overview of pipeline

The pipeline for OPUSeq data analysis is made up of a series of bash scripts. Using Slurm, the master script main.sh submits the other scripts in an order depending on the type of analysis desired. The file var.txt specifies all variables. See the file OPUSeq_pipeline.png for an overview of which scripts will be run depending on specified varibles. The input is paired end FASTQ files which come from either OPUSeq or standard DNA sequencing. The read names must end in *R1.fastq.gz and *R2.fastq.gz. All the .sh and .py scripts, the var.txt file, and the .bed files should be placed in the same subdirectory "scr" inside the main directory. 

###### List of necessary software
- Python (v. 3.7.2)
- Biopython (v. 1.76-py3)
- BWA (v. 0.7.17)
- Bowtie2 (v. 2.3.5.1)
- Samtools (v. 1.14)
- UMI-tools (v. 1.0.0)
- Cutadapt (v. 3.1)
- BamUtil (v. 1.0.15)
- R (v. 4.0.0) 


###### List of variables that need to be specified in the var.txt file
- $maindir: The project directory path.
- $umi_cons: If set to true, consensus sequences will be made using UMIs and boundaries/coordinates. 
- $bound_only_cons: If set to true, consensus sequences will be made using boundaries only. The resulting files will contain "noUMI" in the name. 
Note that $umi_cons and $bound_only_cons are not mutually exclusive. Both types of consensus making can be performed in the same run of the pipeline.
- $dcs: If true, consensus will be made on double strand level. 
- $sscs: If true, consensus will be made on single strand level.
$dcs and $sscs can be set individually, so that e.g. only DCS analysis is performed but not SSCS. If at least one of these is true, then at least one of $umi_cons or $bound_only_cons must be set to true, otherwise an error is raised by main.sh.
- $raw: If true, analyze without making consensus, i.e. just map the fastq and count reference and variant bases per position from the raw BAM files.  
- $library_type: Either "ds_umi" or "no_umi". "ds_umi" is "double strand UMI", i.e. OPUS-seq libraries; "no_umi" is a standard DNA-seq library like KAPA. If $library_type="no_umi", then $umi_cons must be set to false, otherwise an error will be raised by main.sh. 
- $umi_l: Length of UMI (6 for OPUSeq). 
- $adap_err: Maximum error rate when trimming adapters with cutadapt.  
- $min_l: Minimum length of reads left after adapter trimming. If smaller than this, reads are removed. 
- $adap_fw: Path to FASTA file containing forward adapter sequences for trimming with cutadapt. 
- $adap_rv: Path to FASTA file containing reverse adapter sequences for trimming with cutadapt. 
- $edit_dist: The edit distance parameter for UMI tools group. Usually 3. 
- $minmem: Minimal number of family members when making consensus. Usually 3.
- $pileup_q: Minimal quality score for a base to be considered in the mpileup. Usually 30 (not for consensus files, those are not filtered by quality at this stage). 
- $pp_mapq: Minimum mapping quality during filtering with the correct_pair.py script. Usually 36. 
- $bam_trim: Whether to trim BAM files with bamUtil.
- $bam_trim_bp: How many bases to trim from BAM files with bamUtil. 
- $bamutil_path: The path to the folder where bamUtil is installed. 

###### Operations performed by each bash script
- **umi_extract.sh**: Extracts UMIs with UMI tools extract and writes them to the names of each FASTQ read. Resulting files are written to ext_fastq folder. 
- **trim.sh**: Trimming of adapters in two steps. First, the forward adapters are trimmed. Then the reverse (these are only found in molecules where the insert size is small). Trimmed fastq are written to trim_fastq folder (only the final files, the intermediate ones are removed.) Trim log files are also written to this folder.
- **map.sh**: Aligns reads to the human genome (hg19) using BWA aln. The resulting SAM files are sorted, converted to BAM, and indexed using samtools. The number of total, mapped, on target, and flanking reads is calculated using samtools view and the regions specified in HRAS_capture_9.bed and HRAS_capture_9_flank.bed. These stats are written to files *_stats.txt. Output folder is bam. 
- **group.sh**: Corrects UMIs using UMI-tools group (default algorithm). Resulting corrected UMIs are written to TSV files in folder tsv. 
- **pp.sh**: Runs correct_pair.py script on BAM files. Output folder is pp_bam. 
- **consensus.sh**: Runs the script consensus_read.py on the BAM files which have been processed with correct_pair.py. Outputs stats (cmstats folder) and consensus FASTQ files (sscs_fastq and dcs_fastq folders).
- **consensus_map.sh**: Maps the FASTQ files in the dcs_fastq and sscs_fastq folders using bowtie2 in local mode. Output folders are sscs_bam and dcs_bam.
- **bam_filter_trim.sh**: Filters the consensus mapped BAM by mapping quality => 36. Trims all final bam files, if desired (from both ends). The BAM files are rewritten with the new filtered/trimmed ones.
- **pileup.sh**: Produces pileups from files in pp_bam, dcs_bam, sscs_bam using samtools mpileup. The pileup files are written to the pileup folder. Then, the collapse_pileup.py script runs on each pileup file and produces count tables (CSV) with the number of reference and variant bases at each covered position. The pileup files are then removed. 

## Custom Python scripts
###### correct_pair.py
This script takes in BAM files which were produced from paired FASTQ files. If using UMIs is desired, the BAM files should be previously processed by UMI-tools group.

Arguments:
- --input: The path to aligned, paired-end BAM file. 
- --out: Path to output BAM file. 
- --stat: Path to stats file. 
- --UMI-group: Path to .tsv file which is the output from UMI-tools group. 
- --Q: Minimum mapping quality threshold. 

The script filters the aligned paired reads and keeps only those which have the orientation F1R2 or F2R1 (not F1F2 or R1R2) and pass the mapping quality threshold (--Q). If the argument --UMI-group is specified, it also takes the corrected UMI from the BX tag and writes it to the read name. 

###### consensus_read.py 
This Python 3 (v. 3.9.5) script takes in mapped BAM files containing paired reads and creates single-strand consensus sequences (SSCS) and duplex consensus sequences (DCS). These sequences are written out as FASTQ files. In addition, the script produces one txt file with consensus making stats and a few files with tag family stats (optional). The script is based on UnifiedConsensusMaker.py from [https://github.com/Kennedy-Lab-UW/Duplex-Seq-Pipeline](https://github.com/Kennedy-Lab-UW/Duplex-Seq-Pipeline). 

Arguments:
- --input: Path to paired-end BAM file which has been processed with the correct_pair.py script. Required.
- --tagstats: If argument is present, tag family stats will be written out in the form of 1) a tab-separated txt file which contains the number and the fraction of tag families with 1, 2, 3 etc. members; 2) the same values plotted as a barplot in png format; 3) the family sizes of FR:1 vs. RF:2 families from the same DCS. 
- --minmem: Minimum number of reads allowed to comprise a consensus. Default is 3.
- --maxmem: Maximum number of reads allowed to comprise a consensus. Default is 10000. 
- --cutoff: Percentage of nucleotides at a given position in a read that must be identical in order for a consensus (SSCS) to be called at that position. Default is 0.6667.
- --Ncutoff: The maximum fraction of Ns allowed in a consensus sequence. Default is 1, i.e. only sequences with all Ns are filtered out. 
- --write-sscs: If argument is present, print the SSCS reads to file in FASTQ format. 
- --rep_filt: Remove UMIs with homopolymeric runs of nucleotides of this length or longer. Default is 9. 
- --prefix: Sample name to uniquely identify samples. Required.
- --numCores: Number of cores to use for sorting reads. 
- --UMI_in_tag: If argument is present, use UMIs found in the BX tag for consensus making. (UMI-tools group places the corrected UMI in this tag.)
- --UMI_in_name: If argument is present, use UMIs found in the read names for consensus making. Uncorrected UMIs (just after UMI-tools extract) are found in the read names, as well as corrected UMIs after the correct_pair.py script has been run on the file. This is the option used for OPUSeq data analysis.
- --no_output: If argument is present, no SSCS or DCS will be written out (but stats files may be).

UMIs will only be used for consensus making if one of the arguments --UMI_in_tag or --UMI_in_name is provided. 

The script extracts the UMIs from each read in a pair (either from the BX tag or from the read name) and concatenates them into one UMI. It also extracts chromosome number, start coordinate of the forward read, insert size, read pair orientation from the read name (FR or RF, added to read name by the correct_pair.py script), and read number (1 or 2). These are all concatenated and used as the final family tag. (If neither --UMI-in-tag nor --UMI-in-name are specified, then the tag will only contain coordinates and orientation.)

Reads with the same tag are grouped into a family. From each original DNA duplex fragment, four families can be formed: two from the Watson and two from the Crick strand. Watson strand is defined as the one with FR (read 1 forward, read 2 reverse) orientation and so gives rise to FR:1 and FR:2 families. Crick strand has RF orientation (read 1 reverse, read 2 forward) and gives rise to RF:1 and RF:2 families. To ensure that all reads in a family are of the same length, the longest read in a family is found and all other reads are extended with Ns (and their quality strings with 0) until they reach the same length. 

Families with less reads than --minmem or more than --maxmem are excluded. Then, single-strand consensus (SSCS) is made from each remaining family. At each position in the family, bases from all constituent reads are counted. If one base is called in a majority of cases (threshold defined by --cutoff and is set to 0.6667, or 2/3, by default), it is written to the consensus. If no base reaches this threshold, N is written instead. Quality scores written at each position of the consensus sequence are the mean of the Phred quality scores of each base used for consensus. Since Phred scores are logarithmic, it is not mathematically correct to calculate their mean, so this is only a proxy for a quality score. Two FASTQ files with SSCS are written out: one containing the FR:1 and FR:2 consensus sequences and the other - RF:1 and RF:2. The third line in each entry of the FASTQ file contains the family size. 

Further, duplex consensus is called from SSCS pairwise: FR:1 is paired with RF:2 and FR:2 with RF:1. Here, the threshold for consensus is 1 (both bases must be the same), otherwise, N is written. Two FASTQ files are written out, one containing the “read 1” (FR:1 + RF:2) and the other “read 2” (FR:2 + RF:1). The quality score is again calculated as the mean. The third line of each FASTQ entry contains the family sizes of both SSCS families.  

The output file with consensus making statistics contains the number of: 
-	Processed reads
-	Tag families with zero reads
-	Tag families with at least one read, and of these, how many have less than the minimum number of reads (--minmem) and how many have that number or more
-	FR:1 + RF:2 DCS made successfully
-	FR:2 + RF:1 DCS made successfully
-	DCS failed due to missing SSCS
-	DCS excluded due to UMIs with mononucleotide repeats longer than --rep_filt

## collapse_pileup.py
This script takes in a pileup file produced by samtools mpileup from one BAM file and produces a count table in CSV format with counts of reference and variant bases at each position. 

Arguments:
- --input: Path to input pileup file.
- --output: Path to output .csv file. 

The script removes all the special characters in the pileup files before counting, meaning that it ignores indels. It ignores N bases. It also calculates a mean of quality scores of reference and variant bases at each position. This is not mathematically correct due to Phred scores being logarithmic and serves only as a hint about potential differences in quality between reference and variant calls. 

The CSV base count tables were further processed using custom R scripts (see file **OPUSeq.R**).  
