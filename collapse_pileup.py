import csv
import sys
import re
import statistics as stat
from datetime import datetime
from collections import defaultdict
import pandas
from argparse import ArgumentParser

csv.field_size_limit(sys.maxsize)

## This is a dictionary for translating fastq quality score 
## from ASCII to numbers.
qual_trans = {'!':1,'"':2,'#':3,'$':4,'%':5,'&':6,"'":7,'(':8,
              ')':9,'*':10,',':11,'-':12,'.':13,'/':14,'0':15,'1':16,
              '2':17,'3':18,'4':19,'5':20,'6':21,'7':22,'8':23,'9':24,
              ':':25,';':26,'<':27,'=':28,'>':29,'?':30,'@':31,'A':32,
              'B':33,'C':34,'D':35,'E':36,'F':37,'G':38,'H':39,'I':40,
              'J':41,'K':42,'L':43,'M':44,'N':45,'O':46,'P':47,'Q':48,
              'R':49,'S':50,'T':51,'U':52,'V':53,'W':54,'X':55,'Y':56,
              'Z':57,'[':58,'\\':59,']':60,'^':61,'_':62,'`':63,'a':64,
              'b':65,'c':66,'d':67,'e':68,'f':69,'g':70,'h':71,
              'i':72,'j':73,'k':74,'l':75,'m':76,'n':77,'o':78,
              'p':79,'q':80,'r':81,'s':82,'t':83,'u':84,'v':85,
              'w':86,'x':87,'y':88,'z':89,'{':90,'|':91,'}':92,'~':93}
    

def print_time():
    now = datetime.now()
    time = now.strftime("%H:%M:%S")
    print("time =", time)

def qual_calc(base_list,qual_list,qual_dict,base):
## Function to calculate average base quality score for bases 
## aligning to reference or the variant base. 
## Note that this is not an accurate average (since Phred scores are logarithmic) and serves only as a proxy.  
## Provide 'ref' as base to calculate for reference, and 'A','C','G' or 'T' to
## calculate for variants.
    if base == 'ref':
        base1 = '.'
        base2 = ','
    if base == 'A':
        base1 = 'A'
        base2 = 'a'
    if base == 'G':
        base1 = 'G'          
        base2 = 'g'
    if base == 'C':    
        base1 = 'C'
        base2 = 'c'
    if base == 'T':    
        base1 = 'T'
        base2 = 't'
        
    for b in range(0,len(base_list)):
        if base_list[b] == base1:
            qual_dict[base].append(qual_list[b])
        elif base_list[b] == base2:
            qual_dict[base].append(qual_list[b])
    if len(qual_dict[base]) > 0:
        qual_dict[base] = stat.mean(qual_dict[base])
    else:
        qual_dict[base] = 0
    return(qual_dict[base])

def write_df(df_dict,line,base_dict,qual_dict):
## This function makes a dictionary, to be converted to dataframe 
## by pandas, which contains all the info on positions and variants.
## Each key will be a column in the dataframe.
## If a position has no variant, the count is 0 
## and 'NA' is written instead of variant base and quality.
    var_sum = base_dict['A'] + base_dict['C'] + base_dict['G'] + base_dict['T']
    if var_sum == 0:
                df_dict['chr'].append(line[0])
                df_dict['pos'].append(line[1])
                df_dict['ref'].append(line[2])
                df_dict['ref_count'].append(base_dict['ref'])
                df_dict['ref_qual'].append(qual_dict['ref'])
                df_dict['var'].append('NA')
                df_dict['var_count'].append(0)
                df_dict['var_qual'].append('NA')
    if base_dict['A'] != 0:
                df_dict['chr'].append(line[0])
                df_dict['pos'].append(line[1])
                df_dict['ref'].append(line[2])
                df_dict['ref_count'].append(base_dict['ref'])
                df_dict['ref_qual'].append(qual_dict['ref'])
                df_dict['var'].append('A')
                df_dict['var_count'].append(base_dict['A'])
                df_dict['var_qual'].append(qual_dict['A'])
    if base_dict['C'] != 0:
                df_dict['chr'].append(line[0])
                df_dict['pos'].append(line[1])
                df_dict['ref'].append(line[2])
                df_dict['ref_count'].append(base_dict['ref'])
                df_dict['ref_qual'].append(qual_dict['ref'])
                df_dict['var'].append('C')
                df_dict['var_count'].append(base_dict['C'])
                df_dict['var_qual'].append(qual_dict['C'])
    if base_dict['G'] != 0:
                df_dict['chr'].append(line[0])
                df_dict['pos'].append(line[1])
                df_dict['ref'].append(line[2])
                df_dict['ref_count'].append(base_dict['ref'])
                df_dict['ref_qual'].append(qual_dict['ref'])
                df_dict['var'].append('G')
                df_dict['var_count'].append(base_dict['G'])
                df_dict['var_qual'].append(qual_dict['G'])
    if base_dict['T'] != 0:
                df_dict['chr'].append(line[0])
                df_dict['pos'].append(line[1])
                df_dict['ref'].append(line[2])
                df_dict['ref_count'].append(base_dict['ref'])
                df_dict['ref_qual'].append(qual_dict['ref'])
                df_dict['var'].append('T')
                df_dict['var_count'].append(base_dict['T'])
                df_dict['var_qual'].append(qual_dict['T']) 
    return(df_dict)
    
    


def main():
    parser = ArgumentParser()
    parser.add_argument(
        '--input',
        dest = 'in_file',
        required = True,
        help = 'Path to input samtools mpileup file (made from ONE bam file).')
    parser.add_argument(
        '--output',
        dest = 'out_file',
        required = True,
        help = 'Path to output .csv file.')
    arg = parser.parse_args()
    
    counter = 0
    base_dict = {'ref':0,'A':0,'C':0,'G':0,'T':0}
    qual_dict = {'ref':[],'A':[],'C':[],'G':[],'T':[]}
    df_dict = {'chr':[],'pos':[],'ref':[],'ref_count':[],
           'ref_qual':[],'var':[],'var_count':[],'var_qual':[]}
    with open(arg.in_file, 'r', newline = '') as file:                                                             
        file_reader = csv.reader(file, delimiter='\t')
        ##remove all the special characters in the pileup, i.e. those
        ##denoting insertions, deletions, read start and end, etc. 
        ##Leave only characters denoting ref bases, variants and Ns.
        for line in file_reader:
            counter += 1
            line[4] = re.sub(r'\^.|\$|\>|\<', '', line[4])
            ##remove insertions and deletions
            ##of up to 100 bases in length:
            for n in range(0,99):  
                reg_list1 = [r'\+',rf'{n}',r'[A-Za-z]',r'{',rf'{n}',r'}']
                reg1 = re.compile(''.join(reg_list1))
                line[4] = re.sub(reg1, '', line[4])
                reg_list2 = [r'\-',rf'{n}',r'[A-Za-z]',r'{',rf'{n}',r'}']
                reg2 = re.compile(''.join(reg_list2))
                line[4] = re.sub(reg2, '', line[4])
            ##for * character which denotes a deletion already 
            ##mentioned at a previous position, 
            ##we have to delete also the corresponding 
            ##position in the quality string, therefore it is done
            ##this way:
            diff_count = 0
            for match in re.finditer(r'\*+',line[4]):
                star = match.start()
                starend = match.end()
                diff = match.end() - match.start()
                star = star - diff_count
                starend = starend - diff_count
                line[5] = line[5][:star] + line[5][starend:]
                diff_count = diff_count + diff
            line[4] = re.sub(r'\*', '', line[4])
            line[4] = list(line[4])
            qual_list = list(line[5])
            ##translate ASCII to numbers in quality string
            for n in range(0,len(qual_list)):
                qual_list[n] = qual_trans[f'{qual_list[n]}']
            line[5] = qual_list
            ##the following error will come up if some special characters remain in the base list
            if len(line[5]) != len(line[4]):
                print('Error: base and quality lists are of unequal lengths')
                print(line[4])
                print(line[5])
            ##count bases of each kind and calculate average
            ##quality scores for them
            base_dict['ref'] = line[4].count('.')+line[4].count(',')                    
            qual_dict['ref'] = qual_calc(line[4],
                                         line[5],qual_dict,'ref')
            base_dict['A'] = line[4].count('A')+line[4].count('a')                  
            qual_dict['A'] = qual_calc(line[4],
                                       line[5],qual_dict,'A')
            base_dict['C'] = line[4].count('C')+line[4].count('c')                   
            qual_dict['C'] = qual_calc(line[4],
                                       line[5],qual_dict,'C')
            base_dict['G'] = line[4].count('G')+line[4].count('g')                  
            qual_dict['G'] = qual_calc(line[4],
                                       line[5],qual_dict,'G')
            base_dict['T'] = line[4].count('T')+line[4].count('t')                    
            qual_dict['T'] = qual_calc(line[4],
                                       line[5],qual_dict,'T') 
            ##compile the above into a dict
            df_dict = write_df(df_dict,line,base_dict,qual_dict) 
            ##reset qual_dict for next line
            qual_dict = {'ref':[],'A':[],'C':[],'G':[],'T':[]}
    ##convert the dict with all info into a dataframe and write to csv
    final_df = pandas.DataFrame.from_dict(df_dict, orient='columns')
    final_df.to_csv(arg.out_file)

print('Start:')
print_time()

if __name__ == "__main__":
    main()

print('Complete:')
print_time()