import sys
import os
import pysam

class ReadPair(object):
    def __init__(self,read):
        self.id=read.query_name
        self.read1=None
        self.read2=None
        self.f1r2 = 0
        self.r1f2 = 0
        self.f1f2 = 0
        self.r1r2 = 0
        self.update_mate(read)


    def update_mate(self,read):
        if read.is_read1:
            self.read1 = read
        else:
            self.read2 = read

    def update_proper_pair(self, order):
        self.read1.is_proper_pair = True
        self.read2.is_proper_pair = True
        self.id = f"{self.id}#{order}"
        self.read1.query_name = self.read2.query_name = self.id

    def check_proper_pair(self):
        write_output = False
        if self.read1.is_reverse:
            if not self.read2.is_reverse:
                self.r1f2 += 1
                write_output = True
                self.update_proper_pair("RF")
            else:
                self.r1r2 += 1
        else:
            if self.read2.is_reverse:
                self.f1r2 += 1
                write_output = True
                self.update_proper_pair("FR")
            else:
                self.f1f2 += 1
        return write_output, self.f1r2, self.r1f2, self.f1f2, self.r1r2


input_file = pysam.AlignmentFile(sys.argv[1],"rb")
tmp_file = pysam.AlignmentFile("tmp.bam","wb",template = input_file)
read_pairs = dict()
pair_types = [0,0,0,0]

for read in input_file.fetch(until_eof=True):
    if not read.is_unmapped and not read.mate_is_unmapped and read.next_reference_id == read.reference_id:
        name = read.query_name
        if name in read_pairs.keys():
            rp = read_pairs[name]
            rp.update_mate(read)
            write_out, F1R2, R1F2, F1F2, R1R2 = rp.check_proper_pair()
            pair_types = [x + y for x, y in zip(pair_types, [F1R2, R1F2, F1F2, R1R2])]
            if write_out:
                tmp_file.write(rp.read1)
                tmp_file.write(rp.read2)
            else:
                pass
        else:
            read_pairs[name] = ReadPair(read)
input_file.close()
tmp_file.close()

output_bam = sys.argv[1].replace(".bam","_properPair.bam")
pysam.sort("-o",output_bam,"tmp.bam")
pysam.index(output_bam)
os.remove("tmp.bam")

stat_file = open(sys.argv[1].replace(".bam","_readPair.stat"),"w")
stat_file.write("#Read pairs in F1R2: " + str(pair_types[0]) + "\n")
stat_file.write("#Read pairs in R1F2: " + str(pair_types[1]) + "\n")
stat_file.write("#Read pairs in F1F2: " + str(pair_types[2]) + "\n")
stat_file.write("#Read pairs in R1R2: " + str(pair_types[3]) + "\n")
stat_file.close()
