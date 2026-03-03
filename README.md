# Deduper
**A Reference Based PCR Duplicate Removal tool.**

Due to amplification bias, PCR does not always amplify RNA proportionally. This can lead to inaccurate quantification of expression levels of genes/isoforms. Luckily, in many sequencing datasets, UMIs (unique molecular indexes) are attached to cDNA during library preparation prior to amplification, which can help us identify PCR clones. Reads that map to the same genomic location and harbor the same UMI are PCR duplicates, and only one copy must be retained. This tool takes in a SAM file, extracts necessary information from each read (including UMI sequence, chromosome, strand, alignment start site, and soft masking information) to determine which reads are PCR duplicates. A single read out of each set of PCR duplicates (the first one encountered in the file) is retained. All retained reads are written to an output SAM file.

## Running the Tool
Deduper is a Python script ([karumuru_deduper.py](./karumuru_deduper.py)) that can be run from the command line, as follows:
```
./karumuru_deduper.py -f <in.sam> -u <known_umis.txt> -o <out.sam> [-s <dedup_stats.txt>]
```
### Required Inputs
#### `-f` | `--file`
File path to a SAM file that:
- Is sorted (by chromosome, then alignment start position)
    - You may need to use `samtools sort` prior to running this tool 
- Contains each read's UMI sequence at the end of its QNAME, like so:
      ```NS500451:154:HWKTMBGXX:1:11101:15364:1139:GAACAGGT```
- Can contain header lines (starting with `@`)
- Only contains uniquely mapped, high quality reads (the script will not filter out low quality or multi-mapping reads)
- Example: [test.sam](./example_files/test.sam)
#### `-u`  |  `--umi`
File path to a text file that contains list of known UMIs, with each UMI sequence on a new line
- Any reads containing UMIs that aren't in this list will be discarded
- Example: [STL96.txt](./example_files/STL96.txt)
#### `-o`  |  `--outfile`
File path to output deduplicated SAM file
- Will contain all header lines from input SAM file, along with deduplicated reads


### Optional Inputs
##### `s`  |  `--stats`
File path to a output text file containing deduplication statistics
- See [dedup_stats_C1_SE_uniqAlign.txt](./example_files/dedup_stats_C1_SE_uniqAlign.txt) for an example of what this file looks like


