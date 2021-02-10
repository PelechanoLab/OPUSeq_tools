import os
import sys
import pysam
import multiprocessing as mp
import gzip
import datetime
from argparse import ArgumentParser
from collections import defaultdict
from collections import Counter
import statistics
import Bio
from Bio.Seq import Seq


class iteratorWrapper:
    def __init__(self, inIterator, finalValue):
        self.it = inIterator
        self.finalValue = finalValue
        self.endIter = False
    def __iter__(self):
        return(self)
    def __next__(self):
        try:
            temp = next(self.it)
        except StopIteration:
            if self.endIter == False:
                temp = self.finalValue
                self.endIter = True
            else:
                raise(StopIteration)
        return(temp)
    next = __next__

##The consensus caller script takes an object containing all the reads belonging to the same family and makes the consensus sequence. It creates a list containing the counts of each nucleotide at each position, and then for each position, to the consensus seq it writes either the nucleotide that is above a certain threshold/cutoff (0.6667 by default) or N.
    
def consensus_caller(input_reads, cutoff, tag):
    ##first, make all reads the same length by adding Ns until the max length of all reads passed to the function
    lengths = []
    for read in input_reads:
        lengths.append(len(read))
    maxL = max(lengths)
    for read in input_reads:
        if len(read) < maxL:
            diff = maxL - len(read)
            read = read + 'N'*diff
    nuc_identity_list = [0, 0, 0, 0, 0, 0]
    nuc_key_dict = {0: 'T', 1: 'C', 2: 'G', 3: 'A', 4: 'N'}
    consensus_seq = ''

    for i in range(maxL): 
        # Count the types of nucleotides at a position in a read.
        # i is the nucleotide index within a read
        for j in range(len(input_reads)):
        # Do this for every read that comprises a tag family.
        # j is the read index
            try:
                if input_reads[j][i] == 'T':
                    nuc_identity_list[0] += 1
                elif input_reads[j][i] == 'C':
                    nuc_identity_list[1] += 1
                elif input_reads[j][i] == 'G':
                    nuc_identity_list[2] += 1
                elif input_reads[j][i] == 'A':
                    nuc_identity_list[3] += 1
                elif input_reads[j][i] == 'N':
                    nuc_identity_list[4] += 1
                else:
                    nuc_identity_list[4] += 1
                nuc_identity_list[5] += 1
            except Exception:
                break
        ## Calculate which nucleotide at each position passes the threshold; if none pass, write N
        try:
            for j in [0, 1, 2, 3, 4]:
                if (float(nuc_identity_list[j])
                        /float(nuc_identity_list[5])
                        ) >= cutoff:
                    consensus_seq += nuc_key_dict[j]
                    break
                elif j == 4:
                    consensus_seq += 'N'
        except Exception:
            consensus_seq += 'N'
        nuc_identity_list = [0, 0, 0, 0, 0, 0]
        # Reset for the next nucleotide position    
    return consensus_seq

##calculates a quality score from qual_list object as the sum of quality scores of all bases at each position across reads. Should we not calculate some kind of average? Otherwise, won't the score differ a lot based on the family size?
def qual_calc(qual_list):
    #return [round(statistics.mean(qual_score)) for qual_score in zip(*qual_list)]
    Qlengths = []
    for read_qual in qual_list:
        Qlengths.append(len(read_qual))
    maxQL = max(Qlengths)
    for read_qual in qual_list:
        if len(read_qual) < maxQL:
            diffQ = maxQL - len(read_qual)
            read_qual += diffQ * [0]
    return [round(sum(qual_score)/len(qual_score)) for qual_score in zip(*qual_list)]

def main():
    startTime = datetime.datetime.now()
    parser = ArgumentParser()
    parser.add_argument(
        '--input',
        dest = 'in_bam',
        required = True,
        help = 'Path to unaligned, paired-end, bam file.'
        )
    parser.add_argument(
        "--tagstats",
        dest='tagstats',
        action="store_true",
        help="Output tagstats file"
    )
    parser.add_argument(
        '--minmem',
        dest='minmem',
        type=int,
        default=1,
        help="Minimum number of reads allowed to comprise a consensus. [3]"
    )
    parser.add_argument(
        '--maxmem',
        dest='maxmem',
        type=int,
        default=10000,
        help="Maximum number of reads allowed to comprise a consensus. [200]"
    )
    parser.add_argument(
        '--cutoff',
        dest='cutoff',
        type=float,
        default=.6667,
        help=(f"Percentage of nucleotides at a given position "
              f"in a read that must be identical in order "
              f"for a consensus to be called at that position. "
              f"[0.7]"
              )
    )
    parser.add_argument(
        '--Ncutoff',
        dest='Ncutoff',
        type=float,
        default=1,
        help=(f"With --filt 'n', maximum fraction of Ns allowed in a "
              f"consensus [1.0]"
              )
    )
    parser.add_argument(
        '--write-sscs',
        dest='write_sscs',
        action="store_true",
        help="Print the SSCS reads to file in FASTQ format"
    )
    parser.add_argument(
        '--without-dcs',
        dest='without_dcs',
        action="store_true",
        help="Don't print final DCS reads"
    )
    parser.add_argument(
        "--rep_filt",
        action="store",
        type=int,
        dest='rep_filt',
        default=9,
        help=(f"Remove tags with homomeric runs of nucleotides of length "
              f"x. [9]"
              )
    )
    parser.add_argument(
        '--prefix',
        action="store",
        dest='prefix',
        type=str,
        required=True,
        help="Sample name to uniquely identify samples"
    )
    parser.add_argument(
        '--numCores',
        action="store",
        dest="cores",
        type=int,
        default=1,
        help="Number of cores to use for sorting UMI-processed reads."
    )
    parser.add_argument(
        '--numAlignmentReads',
        dest="numAlignReads",
        action="store",
        type=int,
        default=500000,
        help="Number of read pairs to output as fastq to align for determining raw reads on target.  Set to 0 to skip this output.  Will stop when it reaches the end of the file or this number."
    )
    parser.add_argument(
        '--UMI_in_tag',
        dest='UMI_in_tag',
        action="store_true",
        help = "Use UMI for generating consensus sequences."
    )
    parser.add_argument(
        '--UMI_in_name',
        dest='UMI_in_name',
        action="store_true",
        help = "Provided files have UMI in the read name, e.g. OPUS-seq processed with UMI-tools extract."
    )
    parser.add_argument(
        '--no_output',
        dest='no_output',
        action="store_true",
        help = "Do not write DCS bam or fastq files."
    )
    o = parser.parse_args()
    # adjust number of cores
    if o.cores >= mp.cpu_count() and o.cores > 1:
        o.cores = mp.cpu_count() - 1

    in_bam_file = pysam.AlignmentFile(o.in_bam, "rb", check_sq=False)
    temp_bam = pysam.AlignmentFile(f"{o.prefix}.temp.bam",
                                   'wb',template=in_bam_file)

    # Initialize Counters:
    # Counter for reads UMI-processed, and for number of raw reads
    paired_end_count = 0
    # Counter for number of families
    familyCtr = 0
    # Counter for SSCS families where read length is not uniform
    badReadLength = 0
    # Counter for DCS UMIs with bad UMIs
    badUMIs = 0
    # Counter for reads processed
    readsCtr = 0
    # Counter for low familiy size families
    smallFamilySize = 0
    # Counter for unrepresented families
    zeroFamilySize = 0
    # Counter for high-N SSCS filtered
    highN_SSCS = 0
    # Counter for number of SSCS made
    numSSCS = 0
    # counter for number of families that fail to find their partner
    failedDcs = 0
    # Counter for number of DCS made for read 1 (of forward strand)
    numDCS1 = 0
    # Counter for number of DCS made for read 2 (of forward strand)
    numDCS2 = 0
    # Counter for number of high-N DCS filtered
    highN_DCS = 0

    # Open Files
    if o.write_sscs is True:
        if o.UMI_in_tag or o.UMI_in_name:
            sscs_bam_file = pysam.AlignmentFile(f"{o.prefix}.sscs_UMI.bam",
                                                'wb', template=in_bam_file)
        else:
            sscs_bam_file = pysam.AlignmentFile(f"{o.prefix}.sscs.bam",
                                                'wb', template=in_bam_file)
        read1_sscs_fq_file = gzip.open(f"{o.prefix}_read1_sscs.fq.gz", 'wt')
        read2_sscs_fq_file = gzip.open(f"{o.prefix}_read2_sscs.fq.gz", 'wt')

    if o.without_dcs is False:
    #    if o.UMI_in_tag or o.UMI_in_name:
    #        dcs_bam_file = pysam.AlignmentFile(f"{o.prefix}.dcs_UMI.bam",
    #                                           'wb', template=in_bam_file)
    #    else:
    #        dcs_bam_file = pysam.AlignmentFile(f"{o.prefix}.dcs.bam",
    #                                           'wb', template=in_bam_file)
        read1_dcs_fq_file = gzip.open(f"{o.prefix}_read1_dcs.fq.gz", 'wt')
        read2_dcs_fq_file = gzip.open(f"{o.prefix}_read2_dcs.fq.gz", 'wt')

    ## This block of code takes the input file, extracts the coordinates, UMIs and the FR/RF tags, 
    ## and makes a final tag out of all those. If there are no arguments --UMI-in-name or --UMI-in-tag,
    ## only coordinates are used and not UMI. The reads acquire a tag like e.g. "start_insertSize_UMI_FR:1".

    for line in in_bam_file.fetch(until_eof=True):
        original_name = line.query_name
        if line.is_reverse:
            query_name = f"{line.reference_name}_{line.next_reference_start+1}_{abs(line.template_length)}"
        else:
            query_name = f"{line.reference_name}_{line.reference_start+1}_{abs(line.template_length)}"
        pair_tag = line.query_name.split("#")[1]
        if line.is_read1:
            pair_tag = f"{pair_tag}:1"
        else:
            pair_tag = f"{pair_tag}:2"
        ##different ways to handle UMI depending on whether the UMI is in read name (after UMI-tools extract) or in BX tag (after UMI-tools group)
        if o.UMI_in_tag: 
            if o.UMI_in_name:
                raise Exception(f'ERROR: Cannot use UMI from name and tag simultaneously')
            elif "RF" in pair_tag:
                umi_tag = line.get_tag("BX")[6:] + line.get_tag("BX")[0:6]
            else:
                umi_tag = line.get_tag("BX")
            query_name += f"_{umi_tag}#{pair_tag}"            
        elif o.UMI_in_name: 
            label = line.query_name.split("_")[1]
            if "RF" in pair_tag:
                umi_tag = label.split("#")[0][6:] + label.split("#")[0][0:6]
            else:
                umi_tag = label.split("#")[0]
            query_name += f"_{umi_tag}#{pair_tag}"
        else: ##option for data without UMI
            query_name += f"#{pair_tag}"            
        line.query_name = query_name
        line.set_tag('X?', original_name, 'Z') 
        temp_bam.write(line)
        paired_end_count += 1

    in_bam_file.close()
    temp_bam.close()

    pysam.sort("-n", "-@", f"{o.cores}", "-o", f"{o.prefix}.temp.sort.bam", f"{o.prefix}.temp.bam")
    os.remove(f"{o.prefix}.temp.bam")

    seq_dict = {'FR:1': [], 'FR:2': [], 'RF:1': [], 'RF:2': []}
    qual_dict = {'FR:1': [], 'FR:2': [], 'RF:1': [], 'RF:2': []}
    last_seq = {'FR:1': [], 'FR:2': [], 'RF:1': [], 'RF:2': []}
    lengths_dict = {'FR:1': [], 'FR:2': [], 'RF:1': [], 'RF:2': []}
    fam_size_x_axis = []
    fam_size_y_axis = []

    read1_dcs_len = 0
    read2_dcs_len = 0
    in_bam_file = pysam.AlignmentFile(
        f"{o.prefix}.temp.sort.bam", "rb", check_sq=False
    )
    first_line = next(in_bam_file)
    readsCtr += 1
    FinalValue = pysam.AlignedSegment()
    FinalValue.query_name = "FinalValue#FR:1"

    seq_dict[first_line.query_name.split('#')[1]].append(
        first_line.query_sequence
    )
    lengths_dict[first_line.query_name.split('#')[1]].append(
        len(first_line.query_sequence)
    )
    qual_dict[first_line.query_name.split('#')[1]].append(
        list(first_line.query_qualities)
    )
    tag_count_dict = defaultdict(lambda: 0)

    print("Creating consensus reads...")
    
    ########SSCS calling
    for line in iteratorWrapper(in_bam_file.fetch(until_eof=True), FinalValue):
        tag = first_line.query_name.split('#')[0] ##get tag from query name (the UMI and the coordinates)
        ##iterate through the other reads looking for reads with the same tag (coordinate + UMI). When found, write their sequences and quality scores to the dictionaries containing those under the appropriate "subtag".
        if line.query_name.split('#')[0] == tag:
            readsCtr += 1
            subtag = line.query_name.split('#')[1] ##this is either "FR:1", "FR:2", "RF:1", or "RF:2" 
            if subtag == "RF:1" or subtag == "FR:2":
                read_seq = str(Seq(line.query_sequence).reverse_complement())
                qual_list = list(line.query_qualities)
                qual_list.reverse()
            else:
                read_seq = line.query_sequence
                qual_list = list(line.query_qualities)            
            seq_dict[line.query_name.split('#')[1]].append(read_seq)
            qual_dict[line.query_name.split('#')[1]].append(qual_list)
            last_seq[line.query_name.split('#')[1]] = line
        ##count reads with subtags "FR:1", "FR:2", "RF:1", or "RF:2", and if 1&2 numbers don't match, raise error. 
        ##These are reads 1 and 2 of the same orientation, so there should be the same number. 
        else:
            famSizes = {x: len(seq_dict[x]) for x in seq_dict}
            #lengths = []
            if (famSizes['FR:1'] != famSizes['FR:2'] or famSizes['RF:1'] != famSizes['RF:2']):
                raise Exception(f'ERROR: Read counts for Read1 and Read 2 do '
                                f'not match for tag {tag}'
                                )
            for tag_subtype in seq_dict.keys():
                imbalance = False
                #for read in seq_dict[tag_subtype][1:]:
                #    if len(read) != len(seq_dict[tag_subtype][0]):
                #        imbalance = True
                #if imbalance:
                #    l_list = []
                #    for read in seq_dict[tag_subtype]:
                #        l_list.append(len(read))
                #    print(*l_list,sep=" ")
                ##check that all reads within a family are the same length. If not, the family will be empty and counted as badReadLength (see below).
                #for read in seq_dict[tag_subtype]:
                    #lengths.append(len(read))
                #major_length = Counter(lengths).most_common(1)[0][0]
                #for read in seq_dict[tag_subtype]:
                    #if len(read) != major_length:
                        #seq_dict[tag_subtype].remove(read)                        
                if not imbalance:
                    if famSizes[tag_subtype] > 0:
                        tag_count_dict[famSizes[tag_subtype]] += 1 ##dict with counters for each family size, i.e. how many families there are with 1,2,3 etc. members
                        familyCtr += 1
                        ##if family size is zero, count under zeroFamilySize. However, these families are KEPT in the dict. 
                    if famSizes[tag_subtype] == 0:
                        zeroFamilySize += 1
                        ##if family size is smaller than specified minimum (but not zero), delete it from dict and count as smallFamilySize:
                    elif famSizes[tag_subtype] < o.minmem:
                        seq_dict[tag_subtype] = []
                        qual_dict[tag_subtype] = []
                        smallFamilySize += 1
                        ##if family size within the specified limits, run consensus caller on it and write the resulting sequence and family size to the seq_dict; count as numSSCS:
                    elif o.minmem <= famSizes[tag_subtype] <= o.maxmem:
                        # Tag types w/o reads should not be submitted as long as
                        # minmem is > 0
                        seq_dict[tag_subtype] = [
                            consensus_caller(seq_dict[tag_subtype],
                                            o.cutoff,
                                            tag
                                            ),
                            str(famSizes[tag_subtype])
                        ]
                        qual_dict[tag_subtype] = qual_calc(qual_dict[tag_subtype]) ##calculate quality score and write that to the qual_dict
                        numSSCS += 1
                        ##if family size is greater than the specified maximum, call consensus using all the reads up until the limit (the specified maximum). We may want to remove this, or else exclude families with huge sizes since they are likely artifacts.
                    elif famSizes[tag_subtype] > o.maxmem:
                        seq_dict[tag_subtype] = [
                            consensus_caller(seq_dict[tag_subtype][:o.maxmem],
                                            o.cutoff,
                                            tag
                                            ),
                            str(famSizes[tag_subtype])
                        ]
                        qual_dict[tag_subtype] = qual_calc(qual_dict[tag_subtype])
                        numSSCS += 1
                else:
                    badReadLength += 1
                    seq_dict[tag_subtype] = []
                    qual_dict[tag_subtype] = []
            ########Write SSCS to fastq
            if o.write_sscs is True:

                if len(seq_dict['FR:1']) != 0 and len(seq_dict['FR:2']) != 0:
                    FR_read1 = last_seq['FR:1']
                    FR_read2 = last_seq['FR:2']
                    FR_read1.query_sequence = seq_dict['FR:1'][0]
                    FR_read2.query_sequence = seq_dict['FR:2'][0]
                    FR_read1.query_qualities = [x if x < 41 else 41 for x in qual_dict['FR:1']]
                    FR_read2.query_qualities = [x if x < 41 else 41 for x in qual_dict['FR:2']]

                    corrected_qual_score = map(
                        lambda x: x if x < 41 else 41, qual_dict['FR:1']
                    )
                    corrQualStr = ''.join(
                        chr(x + 33) for x in corrected_qual_score
                    )
                    read1_sscs_fq_file.write(f"@{tag}#FR/1 XF:Z:{seq_dict['FR:1'][1]}\n"
                                             f"{seq_dict['FR:1'][0]}\n"
                                             f"+{seq_dict['FR:1'][1]}\n"
                                             f"{corrQualStr}\n"
                                             )

                    corrected_qual_score = map(
                        lambda x: x if x < 41 else 41, qual_dict['FR:2']
                    )
                    corrQualStr = ''.join(
                        chr(x + 33) for x in corrected_qual_score
                    )
                    read2_sscs_fq_file.write(f"@{tag}#FR/2 XF:Z:{seq_dict['FR:2'][1]}\n"
                                             f"{seq_dict['FR:2'][0]}\n"
                                             f"+{seq_dict['FR:2'][1]}\n"
                                             f"{corrQualStr}\n"
                                             )

                    sscs_bam_file.write(FR_read1)
                    sscs_bam_file.write(FR_read2)

                if len(seq_dict['RF:1']) != 0 and len(seq_dict['RF:2']) != 0:

                    RF_read1 = last_seq['RF:1']
                    RF_read2 = last_seq['RF:2']
                    RF_read1.query_sequence = seq_dict['RF:1'][0]
                    RF_read2.query_sequence = seq_dict['RF:2'][0]
                    RF_read1.query_qualities = [x if x < 41 else 41 for x in qual_dict['RF:1']]
                    RF_read2.query_qualities = [x if x < 41 else 41 for x in qual_dict['RF:2']]

                    corrected_qual_score = map(
                        lambda x: x if x < 41 else 41, qual_dict['RF:1']
                    )
                    corrQualStr = ''.join(
                        chr(x + 33) for x in corrected_qual_score
                    )

                    read1_sscs_fq_file.write(f"@{tag}#RF/1 XF:Z:{seq_dict['RF:1'][1]}\n"
                                             f"{seq_dict['RF:1'][0]}\n"
                                             f"+{seq_dict['RF:1'][1]}\n"
                                             f"{corrQualStr}\n"
                                             )

                    corrected_qual_score = map(
                        lambda x: x if x < 41 else 41, qual_dict['RF:2']
                    )
                    corrQualStr = ''.join(
                        chr(x + 33) for x in corrected_qual_score
                    )

                    read2_sscs_fq_file.write(f"@{tag}#RF/2 XF:Z:{seq_dict['RF:2'][1]}\n"
                                             f"{seq_dict['RF:2'][0]}\n"
                                             f"+{seq_dict['RF:2'][1]}\n"
                                             f"{corrQualStr}\n"
                                             )

                    sscs_bam_file.write(RF_read1)
                    sscs_bam_file.write(RF_read2)
                    
    
            if o.without_dcs is False: ########DCS calling
                if len(seq_dict['FR:1']) != 0 and len(seq_dict['RF:2']) != 0:
                    numDCS1 += 1
                    ##dcs_read is a list of the DCS made by consensus caller and the family sizes of SSCSs
                    dcs_read_1 = [
                        consensus_caller(
                            [seq_dict['FR:1'][0], seq_dict['RF:2'][0]],
                            1, ##"1" means that a base is called ONLY if it is the same in both the SSCSs, otherwise N is written. 
                            tag,
                        ),
                        seq_dict['FR:1'][1], seq_dict['RF:2'][1]
                    ]
                    dcs_read_1_qual = qual_calc([qual_dict['FR:1'], qual_dict['RF:2']])
                    read1_dcs_len = len(dcs_read_1[0])
                    fam_size_x_axis.append(int(seq_dict['FR:1'][1])) ##this is for plotting family sizes later
                    fam_size_y_axis.append(int(seq_dict['RF:2'][1]))
                    ##check that there are not more Ns that acceptable (the N cutoff specified). If there are too many, re-write everything with N and quality with zero.
					##The default cutoff is 1, which means that this filter is NOT applied at all by default.
                    if dcs_read_1[0].count('N') / read1_dcs_len > o.Ncutoff:
                        highN_DCS += 1
                        dcs_read_1[0] = 'N' * read1_dcs_len
                        dcs_read_1_qual = [0 for x in range(read1_dcs_len)]
                else:
                    failedDcs += 1 ##if one or both of the strands are not present
                if len(seq_dict['RF:1']) != 0 and len(seq_dict['FR:2']) != 0:
                    numDCS2 += 1
                    dcs_read_2 = [
                        consensus_caller(
                            [seq_dict['RF:1'][0], seq_dict['FR:2'][0]],
                            1,
                            tag,
                        ),
                        seq_dict['RF:1'][1], seq_dict['FR:2'][1]
                    ]
                    dcs_read_2_qual = qual_calc([qual_dict['FR:2'], qual_dict['RF:1']])
                    read2_dcs_len = len(dcs_read_2[0])
                    if dcs_read_2[0].count('N') / read2_dcs_len > o.Ncutoff:
                        highN_DCS += 1
                        dcs_read_2[0] = 'N' * read2_dcs_len
                        dcs_read_2_qual = [0 for x in range(read2_dcs_len)]
                else:
                    failedDcs += 1
                ##check for homopolymer repeats in UMI
                if (read1_dcs_len != 0
                        and read2_dcs_len != 0
                        and tag.count('N') == 0
                        and 'A' * o.rep_filt not in tag
                        and 'C' * o.rep_filt not in tag
                        and 'G' * o.rep_filt not in tag
                        and 'T' * o.rep_filt not in tag
                ):
                ##write DCS fastq files
                    if not o.no_output:
                        r1QualStr = ''.join(chr(x + 33) for x in dcs_read_1_qual)
                        r2QualStr = ''.join(chr(x + 33) for x in dcs_read_2_qual)
                        read1_dcs_fq_file.write(
                            f"@{tag}/1 XF:Z:{dcs_read_1[1]}:{dcs_read_1[2]}\n"
                            f"{dcs_read_1[0]}\n"
                            f"+{dcs_read_1[1]}:{dcs_read_1[2]}\n" ##this writes the family sizes of the SSCSs to the third line in fastq files
                            f"{r1QualStr}\n"
                        )
                        read2_dcs_fq_file.write(
                            f"@{tag}/2 XF:Z:{dcs_read_2[1]}:{dcs_read_2[2]}\n"
                            f"{dcs_read_2[0]}\n"
                            f"+{dcs_read_2[1]}:{dcs_read_2[2]}\n"
                            f"{r2QualStr}\n"
                        )
                    ##write DCS bam file
                        #read1 = last_seq['FR:1']
                        #read2 = last_seq['FR:2']
                        #read1.query_name = tag
                        #read2.query_name = tag
                        #read1.query_sequence = dcs_read_1[0]
                        #read2.query_sequence = dcs_read_2[0]
                        #read1.query_qualities = [x for x in dcs_read_1_qual]
                        #read2.query_qualities = [x for x in dcs_read_2_qual]
                        #dcs_bam_file.write(read1)
                        #dcs_bam_file.write(read2)


                elif (tag.count('N') != 0
                      or 'A' * o.rep_filt in tag
                      or 'C' * o.rep_filt in tag
                      or 'G' * o.rep_filt in tag
                      or 'T' * o.rep_filt in tag
                ):
                    badUMIs += 1 ##badUMIs counts UMIs with either homopolymer repeats OR an N.
                    #print(tag)

            if line != FinalValue:
                readsCtr += 1
                # reset conditions for next tag family
                first_line = line
                seq_dict = {'FR:1': [], 'FR:2': [], 'RF:1': [], 'RF:2': []}
                qual_dict = {'FR:1': [], 'FR:2': [], 'RF:1': [], 'RF:2': []}
                read1_dcs_len = 0
                read2_dcs_len = 0
                dcs_read_1 = ''
                dcs_read_2 = ''

                seq_dict[line.query_name.split('#')[1]].append(line.query_sequence)
                # Now add initializing data for new tag
                qual_dict[first_line.query_name.split('#')[1]].append(
                    list(first_line.query_qualities)
                )

# Try to plot the tag family sizes
    if o.tagstats is True:
        tag_stats_file = open(o.prefix + ".tagstats.txt", 'w')

        x_value = []
        y_value = []
        total_reads = sum(
            [tag_count_dict[tag_family_size]
             * tag_family_size for tag_family_size in tag_count_dict.keys()
             ])

        for tag_family_size in sorted(tag_count_dict.keys()):
            fraction = (tag_count_dict[tag_family_size] * tag_family_size)
            fraction /= float(total_reads)
            tag_stats_file.write(
                f'{tag_family_size}\t'
                f'{tag_count_dict[tag_family_size]}\t'
                f'{fraction}\n'
                )
            x_value.append(tag_family_size)
            y_value.append(fraction)

        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            plt.figure(1)
            plt.bar(x_value, y_value)
            plt.xlabel('Family Size')
            plt.ylabel('Proportion of Total Reads')
            plt.savefig(f"{o.prefix}_family_size.png",
                        bbox_inches='tight'
                        )
            plt.figure(2)
            if len(fam_size_x_axis) != 0:
                plt.scatter(fam_size_x_axis, fam_size_y_axis, alpha=.1)
                plt.xlabel('Family size for FR:1')
                plt.ylabel('Family size for RF:2')
                plt.xlim(0, max(fam_size_x_axis))
                plt.ylim(0, max(fam_size_y_axis))
            plt.savefig(f"{o.prefix}_fam_size_relation.png",
                        bbox_inches='tight'
                        )

        except ImportError:
            sys.stderr.write(
                'matplotlib not present. Only tagstats file will be generated.'
                )

        tag_stats_file.close()
    endTime = datetime.datetime.now()
    # Print consensus making statistics
    startTimeStr = startTime.strftime("%A, %d. %B %Y %I:%M%p")
    endTimeStr = endTime.strftime("%A, %d. %B %Y %I:%M%p")
    cmStatsFile = open(f"{o.prefix}_cmStats.txt", 'w')
    cmStatsFile.write(
        f"Consensus Making Statistics:\n"
        f"Command: {' '.join(sys.argv)}\n"
        f"Started at {startTimeStr}\n"
        f"Finished at {endTimeStr}\n"
        f"{paired_end_count} reads UMI processed\n"
        f"{readsCtr} reads processed\n"
        f"{badReadLength} families with non-uniform read length, excluded\n"
        f"{zeroFamilySize} families with zero reads\n"
        f"{familyCtr} families with uniform reads, of these:\n"
        f"\t{smallFamilySize} families with family size < {o.minmem}\n"
        f"\t{numSSCS} families with size > {o.minmem}, made into SSCS\n"
        f"{numDCS1} FR:1+RF:2 DCS made successsfully\n"
        f"{numDCS2} FR:2+RF:1 DCS made successsfully\n"
        f"{failedDcs} DCS failed due to missing SSCS\n"
        f"{badUMIs} DCS excluded due to UMIs with mononucleotide repeats or Ns\n"
        )
    cmStatsFile.close()
    sys.stderr.write(
        f"Consensus Making Statistics:\n"
        f"Command: {' '.join(sys.argv)}\n"
        f"Started at {startTimeStr}\n"
        f"Finished at {endTimeStr}\n"
        f"{paired_end_count} reads UMI processed\n"
        f"{readsCtr} reads processed\n"
        f"{badReadLength} families with non-uniform read length, excluded\n"
        f"{zeroFamilySize} families with zero reads\n"
        f"{familyCtr} families with uniform reads, of these:\n"
        f"\t{smallFamilySize} families with family size < {o.minmem}\n"
        f"\t{numSSCS} families with size > {o.minmem}, made into SSCS\n"
        f"{numDCS1} FR:1+RF:2 DCS made successsfully\n"
        f"{numDCS2} FR:2+RF:1 DCS made successsfully\n"
        f"{failedDcs} DCS failed due to missing SSCS\n"
        f"{badUMIs} DCS excluded due to UMIs with mononucleotide repeats or Ns\n"
        )
    ##remove the temp sorted bam    
    os.remove(f"{o.prefix}.temp.sort.bam")    



if __name__ == "__main__":
    main()