#!/usr/bin/env python
import argparse
import re

# set global variables to hold inputs
def get_args():
    parser = argparse.ArgumentParser(description="A script that deduplicates a SAM file by removing PCR duplicates.")
    parser.add_argument("-f", help="designates absolute file path to sorted SAM file", required=True, type=str)
    parser.add_argument("-o", help="designates absolute file path to deduplicated SAM file", required=True, type=str)
    parser.add_argument("-u", help="designates file containing list of UMIs", required=True, type=str)
    return parser.parse_args()
args = get_args()

# create a set containing known UMIs
umis = set()
with open(args.u, 'r') as u_fh:
    for line in u_fh:
        umis.add(line.strip('\n'))

with open(args.f, 'r') as i_fh, open(args.o, 'w') as o_fh:
    combo_set = set()
    current_chrom = None
    for line in i_fh:
        # write line to output file if it is a header
        if line[0] == "@":
            o_fh.write(line)
            pass 
        
        # split line into list to access elements
        split_line = line.strip("\n").split("\t")

        # check if UMI is a known UMI
        umi = split_line[0][-8:]
        if umi in umis:

            # empty combo set if all reads for a chromosome have been parsed through
            chrom = split_line[2]
            if chrom != current_chrom:
                combo_set = set()
            
            # determine strand of read
            bf = int(split_line[1])
            if ((bf & 16) == 16):
                strand = "-"
            else:
                strand = "+"

            # determine 5' start position of read
            map_start_pos = int(split_line[3])
            cigar = line[5]
            
            ## plus strand calculation
            if strand == "+":
                soft_clip = re.match("\d+S", cigar)
                if soft_clip:
                    clip_amount = int(soft_clip.group()[:-1])
                    start_pos_5p = map_start_pos - clip_amount
                else:
                    start_pos_5p = map_start_pos
            
            ## minus strand calculation
            else:
                cigar_split = re.findall(r'\d+[MNDS]', cigar)
                # ignore soft clipping at beginning
                if 'S' in cigar_split[0]:
                    cigar_nums = [int(chars[:-1]) for chars in cigar_split[1:]]
                else:
                    cigar_nums = [int(chars[:-1]) for chars in cigar_split]
                cigar_sum = sum(cigar_nums)
                start_pos_5p = map_start_pos + cigar_sum

            # check if UMI-strand-5' start position combination has been encountered before
            combo = f"{umi}-{strand}-{start_pos_5p}"
            if combo not in combo_set:
                o_fh.write(line)
                combo_set.add(combo)




            
                
            


