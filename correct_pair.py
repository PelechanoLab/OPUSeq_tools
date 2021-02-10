import sys
import os
import pysam
import csv
from argparse import ArgumentParser

class ReadPair(object):
    def __init__(self,read):
        self.id=read.query_name
        self.read1=None ##this is updated with the read once update_mate runs
        self.read2=None
        self.f1r2 = 0
        self.r1f2 = 0
        self.f1f2 = 0
        self.r1r2 = 0
        self.BX = None
        self.update_mate(read)

    def update_UMI_group(self, umi_dict):
        self.BX = umi_dict[self.id]

    def update_mate(self,read):
        if read.is_read1:
            self.read1 = read
        else:
            self.read2 = read

    ##This function 1) flags the reads as a proper pair, which may not have been done previously by the aligner; 2) adds a tag to the end of the read name/ID, either "FR" if reads are F1R2 or "RF" if reads are F2R1.         
    def update_proper_pair(self, order): 
        self.read1.is_proper_pair = True
        self.read2.is_proper_pair = True
        self.id = f"{self.id}#{order}"
        self.read1.query_name = self.read2.query_name = self.id
        self.read1.set_tag("BX", self.BX)
        self.read2.set_tag("BX", self.BX)

    def check_proper_pair(self):
        write_output = False
        if self.read1.is_reverse:
            if not self.read2.is_reverse: ##If read 1 is reverse and read 2 is forward, count the pair as F2R1, tag it, and mark for inclusion in the output bam.
                self.r1f2 += 1
                write_output = True
                self.update_proper_pair("RF")
            else:
                self.r1r2 += 1 ##if read 1 is reverse and read 2 is also reverse, count it as R1R2 and DO NOT mark for inclusion in final bam.
        else:
            if self.read2.is_reverse: ##If read 2 is reverse and read 1 is forward, count pair as F1R2, tag it, and mark for inclusion.
                self.f1r2 += 1 
                write_output = True
                self.update_proper_pair("FR")
            else:
                self.f1f2 += 1 ##if read 1 is forward and read 2 is forward, mark as F1F2 and DO NOT mark for inclusion.
        ##write out the pairs which are either F1R2 or F2R1 to the output; write out the counters for classifying reads
        return write_output, self.f1r2, self.r1f2, self.f1f2, self.r1r2

def main():
    parser = ArgumentParser()
    parser.add_argument(
        '--input',
        dest = 'in_bam',
        required = True,
        help = 'Path to aligned, paired-end, bam file.'
        )
    parser.add_argument(
        '--out',
        dest = 'out_bam',
        required = True,
        help = 'Path to output bam file.'
        )
    parser.add_argument(
        '--stat',
        dest = 'stats_file',
        required = True,
        help = 'Path to stats file.'
        )
    parser.add_argument(
        "--UMI-group",
        dest='umi_group',
        help="Path to grouped UMI file (tsv format)."
    )

    arg = parser.parse_args()
    input_file = pysam.AlignmentFile(arg.in_bam,"rb")
    tmp_file = pysam.AlignmentFile("tmp.bam","wb",template = input_file)
    read_pairs = dict()
    umi_group_dict = dict()
    pair_types = [0,0,0,0]

    if arg.umi_group:
        group_file = open(arg.umi_group)
        group_info = csv.reader(group_file, delimiter="\t")
        umi_group_dict = {x[0]:x[6] for x in group_info}
        group_file.close()

    ##When a read name is first encountered and gets added to read_pairs dict, does the loop then go over it again? Why? If not, how do both read 1 and read 2 get added to the ReadPairs object?
    for read in input_file.fetch(until_eof=True):
        ##if both read and mate are mapped and their IDs match
        if not read.is_unmapped and not read.mate_is_unmapped and read.next_reference_id == read.reference_id:
            name = read.query_name
            if name in read_pairs.keys(): ##if read name already in dict
                rp = read_pairs[name] ##rp is now the ReadPairs class object
                ##writes the read to either read1 or read2 attribute of the ReadPairs object, depending on whether it is read 1 or 2.
                rp.update_mate(read)
                ##check_proper_pair will classify reads as F1R2, F2R1, F1F2, or R1R2 and mark for inclusion (write_out = TRUE) only the former two categories. It will also have added an FR or RF tag to the end of the read name.
                write_out, F1R2, R1F2, F1F2, R1R2 = rp.check_proper_pair()
                pair_types = [x + y for x, y in zip(pair_types, [F1R2, R1F2, F1F2, R1R2])] ##pair_types is a counter for the type of pair, for each read, one of these four categories will be a 1 after check_proper_pair and the others will be 0, so this 1 is hereby added to the proper position in counter.
                if write_out: ##write out the reads marked for it to bam.
                    tmp_file.write(rp.read1)
                    tmp_file.write(rp.read2)
                else:
                    pass
            else: ##if read is not in dict yet, add it to read_pairs dict by converting it into a ReadPair class object.
                newRP = ReadPair(read)
                if arg.umi_group:
                    newRP.update_UMI_group(umi_group_dict)
                read_pairs[name] = newRP
    input_file.close()
    tmp_file.close()

    output_bam = arg.out_bam
    pysam.sort("-o",output_bam,"tmp.bam")
    pysam.index(output_bam)
    os.remove("tmp.bam")

    stat_file = open(arg.stats_file,"w")
    stat_file.write("#Read pairs in F1R2: " + str(pair_types[0]) + "\n")
    stat_file.write("#Read pairs in R1F2: " + str(pair_types[1]) + "\n")
    stat_file.write("#Read pairs in F1F2: " + str(pair_types[2]) + "\n")
    stat_file.write("#Read pairs in R1R2: " + str(pair_types[3]) + "\n")
    stat_file.close()

if __name__ == "__main__":
    main()