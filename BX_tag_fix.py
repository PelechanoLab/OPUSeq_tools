import sys
import os
import pysam


class ReadPair(object):
    def __init__(self,read):
        self.id=read.query_name
        self.read1=None ##this is updated with the read once update_mate runs
        self.read2=None
        self.BX_tag=None
        self.update_mate(read)

    def update_mate(self,read):
        if read.is_read1:
            self.read1 = read
            self.BX_tag = read.get_tag("BX")
        else:
            self.read2 = read
            
    def set_tag(self,read):
        read.set_tag("BX",self.BX_tag)


input_file = pysam.AlignmentFile(sys.argv[1],"rb")
tmp_file = pysam.AlignmentFile("tmp.bam","wb",template = input_file)

read_pairs = dict()

for read in input_file.fetch(until_eof=True):
    name = read.query_name
    if name not in read_pairs.keys():
        read_pairs[name] = ReadPair(read)
    else:
        rp = read_pairs[name]
        rp.update_mate(read)
        #print(rp.BX_tag)
        if read.is_read2:
            rp.set_tag(read)
            #print(read.get_tag("BX"))
        tmp_file.write(rp.read1)
        tmp_file.write(rp.read2)

input_file.close()
tmp_file.close()

output_bam = sys.argv[1].replace(".bam","_tagfix.bam")
pysam.sort("-o",output_bam,"tmp.bam")
pysam.index(output_bam)
os.remove("tmp.bam")