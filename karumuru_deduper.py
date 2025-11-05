#!/usr/bin/env python
import argparse
import re


### FUNCTIONS ----------------------------------------------------------------------------

# set global variables to hold inputs
def get_args():
    parser = argparse.ArgumentParser(description="A script that deduplicates a SAM file by removing PCR duplicates.")
    parser.add_argument("-f", help = "designates absolute file path to sorted SAM file", required = True, type = str)
    parser.add_argument("-o", help = "designates absolute file path to output deduplicated SAM file", required = True, type = str)
    parser.add_argument("-u", help = "designates file containing list of UMIs", required = True, type = str)
    parser.add_argument("-s", help = "[OPTIONAL] designates absolute file path to output deduplication statistics file", required = False, type = str)
    return parser.parse_args()
args = get_args()


# define a function to determine which strand a read is located on
def determine_strand(SAM_line):
    """
    A function that determines which strand a read is located on.

    INPUT (list): A SAM line split by tabs into a list, with order of fields in SAM file retained and tab and new line characters removed. Each item in list is a string.
    OUTPUT (string): "+" (denoting forward strand) or "-" (denoting minus strand)

    EXAMPLE INPUT: ['NS500451:154:HWKTMBGXX:2:21202:19092:19786:ACACTGTG', '16', '1', '15911793', '255', '72M', '*', '0', '0', 'ATTTGGGGCTATGTGAGTAAGGTGCTCAAAGTGTTGGGTTCCTGGGAGCTGGAGTTACAGGCATTTGTGAGC', 'EEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEA<EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/', 'NH:i:1', 'HI:i:1', 'AS:i:70', 'nM:i:0']
    """

    # extract bitwise flag from line
    bf = int(SAM_line[1])
    if ((bf & 16) == 16):
        strand = "-"
    else:
        strand = "+"

    return(strand)

# define a function to calculate 5' start position of a read
def determine_start_pos(SAM_line, strand):
    """
    A function that calculates 5' start position of a read by parsing alignment information from the read's CIGAR string. Soft clipping is taken into account, and calculation is strand-dependent.
    
    INPUT 1 (list): A SAM line split by tabs into a list, with order of fields in SAM file retained and tab and new line characters removed. Each item in list is a string.
    INPUT 2 (string): Strandedness of read. Either "+" (denoting forward strand) or "-" (denoting minus strand).
    OUTPUT (integer): The 5' start position of read (1 based). Note that chromosome is NOT outputted - only position along chromosome.

    For reads on the forward strand, any soft clipping at the beginning of the read is subtracted from the location that the read begins mapping to the reference.
    For reads on the reverse strand, the length of the reference that is consumed by the read is added to the location that the read begins mapping to the reference. Soft clipping at the beginning of the read is disregarded.

    EXAMPLE INPUT 1: ['NS500451:154:HWKTMBGXX:2:21202:19092:19786:ACACTGTG', '16', '1', '15911793', '255', '72M', '*', '0', '0', 'ATTTGGGGCTATGTGAGTAAGGTGCTCAAAGTGTTGGGTTCCTGGGAGCTGGAGTTACAGGCATTTGTGAGC', 'EEEEEEEEEEEEEEEEEEEEEEAEEEEEEEEEEEEEEEA<EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE/', 'NH:i:1', 'HI:i:1', 'AS:i:70', 'nM:i:0']
    """

    map_start_pos = int(SAM_line[3])
    cigar = SAM_line[5]
    
    # plus strand calculation
    if strand == "+":
        # look for soft clipping at begommomg of cigar string
        soft_clip = re.match(r'\d+S', cigar)
        if soft_clip:
            clip_amount = int(soft_clip.group()[:-1])
            # subtract soft clip amount from 5' mapping start position
            start_pos_5p = map_start_pos - clip_amount
        else:
            # no adjustment necessary if there is no soft clipping at beginning of read
            start_pos_5p = map_start_pos
    
    # minus strand calculation
    if strand == "-":
        # get a list of number-letter duos that consume the reference in CIGAR string
        cigar_split = re.findall(r'\d+[MNDS]', cigar)
        
        # ignore soft clipping at beginning
        if 'S' in cigar_split[0]:
            # extract numbers from CIGAR string
            cigar_nums = [int(chars[:-1]) for chars in cigar_split[1:]]
        else:
            # extract numbers from CIGAR string
            cigar_nums = [int(chars[:-1]) for chars in cigar_split]
        # get length of reference consumed by read
        cigar_sum = sum(cigar_nums)
        # add length to 5' mapping start position to get 5' start position
        start_pos_5p = map_start_pos + cigar_sum

    return start_pos_5p



### INITIAL SET UP ----------------------------------------------------------------------------

# create a set containing known UMIs
umis = set()
with open(args.u, 'r') as u_fh:
    for line in u_fh:
        umis.add(line.strip('\n'))

# initialize a dictionary for counting unique reads per chromosome
chr_reads = {}

# initialize counters for deduplication summary statistics
num_reads = 0
headers = 0
unknown_umis = 0
unique_reads = 0
dup_reads = 0



### PARSING SAM FILE ----------------------------------------------------------------------------

# open input SAM file for reading and output SAM file for writing
with open(args.f, 'r') as i_fh, open(args.o, 'w') as o_fh:

    # initialize a tracker for current chromosome number/name
    current_chrom = None

    # loop through lines of input SAM file
    for line in i_fh:
        # write line to output file if it is a header
        if line[0] == "@":
            o_fh.write(line)
            headers += 1
            # immediately skip to next line in input SAM file
            continue
        
        num_reads += 1

        # split line into list to access elements
        split_line = line.strip("\n").split("\t")

        # check if UMI is a known UMI
        umi = split_line[0][-8:]
        if umi in umis:
            # extract chromosome name/number
            chrom = split_line[2]
            # initialize variables when a new chromosome is encountered
            if chrom != current_chrom:
                # initalize a set to contain UMI-strand-5' start position combinations
                combo_set = set()
                chr_reads[chrom] = 0
                current_chrom = chrom
                
            # determine strand of read
            strand = determine_strand(split_line)

            # determine 5' start position of read
            start_pos_5p = determine_start_pos(split_line, strand)

            # concatenate read information into a single string
            combo = f"{umi}-{strand}-{start_pos_5p}"
            # check if UMI-strand-5' start position combination has been encountered before
            if combo not in combo_set:
                unique_reads += 1
                # write read to output SAM file
                o_fh.write(line)
                combo_set.add(combo)
                chr_reads[chrom] += 1
            else:
                # don't write read to output SAM file if it is a duplicate
                dup_reads += 1

        else:
            unknown_umis += 1



# SUMMARIZING RESULTS ----------------------------------------------------------------------------

if args.s:

    # define a function to sort chromosome names in a human readable order
    def sorted_nicely(ls): 
        """ 
        Sort a list numerically, then alphabetically. 
        All numeric values are listed first, sorted numerically. Then all alphabet/alphanumeric values are listed, sorted alphabetically.

        INPUT (list): a list of strings (can be a mix of numeric, alphabet, and alphanumeric strings)
        OUTPUT (list): a sorted list of strings
        """ 
        convert_type = lambda text: int(text) if text.isdigit() else text 
        alphanum_key = lambda key: [convert_type(c) for c in re.split('([0-9]+)', key)] 
        return sorted(ls, key = alphanum_key)
    
    sorted_chrom_names = sorted_nicely(list(chr_reads.keys()))

    # write summary statistics to output file
    with open(args.s, 'w') as s_fh:
        s_fh.write("Deduplication Summary Statistics\n\n")
        s_fh.write(f"Total Input Reads: {num_reads}\n")
        s_fh.write(f"Header Lines: {headers}\n")
        s_fh.write(f"Unique Reads: {unique_reads}\n")
        s_fh.write(f"Duplicate (Removed) Reads: {dup_reads}\n")
        s_fh.write(f"Reads with Unnkown UMIs: {unknown_umis}\n")
        s_fh.write(f"\nUnique Reads per Chromosome:\n")
        # list numerical chromsomes first, followed by alphanumeric/alpha chromosomes
        for key in sorted_chrom_names:
            s_fh.write(f"{key}\t{chr_reads[key]}\n")


