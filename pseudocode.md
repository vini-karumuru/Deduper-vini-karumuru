### Define the Problem

Due to amplification bias, PCR does not always amplify RNA proportionally. This can lead to inaccurate quantification of expression levels of genes/isoforms. Luckily, in the sequencing dataset we are given, a UMI (unique molecular index) was attached to each piece of cDNA during library preparation prior to amplification, and this can help us identify PCR clones. Reads that map to the same genomic location and harbor the same UMI are PCR duplicates, and only one copy must be retained. We are given a SAM file, which we must parse through to extract necessary information from each read (including UMI sequence, chromosome, strand, alignment start site, soft masking information, etc.) to determine which reads are PCR duplicates. Then we must retain a single copy out of each set of PCR duplicates (the first one encountered in the file). All retained reads are to be written to an output SAM file.

### Pseudocode

- Use argparse for input variables to script:
    - `-f`, `—file`: designates absolute file path to sorted sam file
    - `-o`, `—outfile`: designates absolute file path to deduplicated sam file
    - `-u`, `—umi`: designates file containing list of UMIs
    - `-h`, `—help`: prints a help message explaining options
- open UMI file for reading:
    - initialize an empty list that will contain UMI sequences
    - loop through each line + append UMI sequence to list
    - turn list into a set (for faster access in future) called umis
- open input SAM file for reading:
    - open output SAM file for writing:
        - initialize a current_location variable to be an empty string
        - initialize a current_reads_plus variable to be an empty list
        - initialize a current_reads_minus variable to be an empty list
        - use a for loop to loop through each line of input SAM file
            - strip line of new line character
            - run check_umi function on the line → if output is “invalid” :
                - continue in loop so that the next line will be read in
            - run check_location function on the line → if output of check_location is equal to current_location:
                - run check_strandedness function on the line → if output is equal to “plus”, add to current_reads_plus, else add to current_reads_minus
            - if not:
                - all reads for current location have been collected, so store the line and location in next_read and next_location variable
                - initialize a variable called current_umis to be an empty list
                - loop through current_reads_plus:
                    - check UMI → if it is already existing current_umis, do nothing, else add it to current_umis and write read to output file
                - empty current_umis list
                - repeat for current_reads_minus
                - overwrite current_location variable to be value of next_location
                - empty current_reads_plus
                - empty current_reads_minus
                - check_strand of next_read, and add it to current_reds_plus or current_reads_minus depending on output

### High Level Functions

```bash
check_umi():
		'''This function takes in a read as a string and extracts the UMI, checks it against list of known UMIs, then returns "valid" if it is known and "invalid" if it is not'''
		return(validity)
		
# Example input:
check_umi("NS500451:154:HWKTMBGXX:1:11101:1232:1273:GGATAACG	0	2	52308305	36	71M	*	0	0	CTTCGTATTCCTGGGTAATGGTCTGGGGGAAGAAGCCTTTGCCTCTGTCCTCTTCATACTCAGCTTTGTAG	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<<EEEEEEEEEEA<EEEEEEEEE<EEEEAAEE<EEEAEA	MD:Z:4A66	NH:i:1	HI:i:1	NM:i:1	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU")

# Expected output:
"valid"
```

```bash
check_location():
		'''This function takes in a read as a string and parses necessary information (chromosome, starting alignment location, soft masking amount, strand) and calculates and outputs chromosome and genomic start location.'''
		return(location)
		
# Example input:
check_location("NS500451:154:HWKTMBGXX:1:11101:1232:1273:GGATAACG	0	2	52308305	36	71M	*	0	0	CTTCGTATTCCTGGGTAATGGTCTGGGGGAAGAAGCCTTTGCCTCTGTCCTCTTCATACTCAGCTTTGTAG	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<<EEEEEEEEEEA<EEEEEEEEE<EEEEAAEE<EEEAEA	MD:Z:4A66	NH:i:1	HI:i:1	NM:i:1	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU")

# Expected output:
"2-52308305"

```

```bash
check_strand():
		'''This function takes in a read as a string and checks the bitwise flag and returns "plus" or "minus"'''
		return(strand)
		
# Example input:
check_strand("NS500451:154:HWKTMBGXX:1:11101:1232:1273:GGATAACG	0	2	52308305	36	71M	*	0	0	CTTCGTATTCCTGGGTAATGGTCTGGGGGAAGAAGCCTTTGCCTCTGTCCTCTTCATACTCAGCTTTGTAG	6/EEEEEEEEEEEEEEEEEEEEEEEEEEEEEE<<EEEEEEEEEEA<EEEEEEEEE<EEEEAAEE<EEEAEA	MD:Z:4A66	NH:i:1	HI:i:1	NM:i:1	SM:i:36	XQ:i:40	X2:i:0	XO:Z:UU")

# Expected output:
"plus"
```