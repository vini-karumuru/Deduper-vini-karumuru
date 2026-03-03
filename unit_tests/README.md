# Deduper Unit Tests

## File Descriptions
- `test_input.sam` is the original, unsorted SAM file containing test cases for Deduper
- `test_input_sorted.sam` is the sorted version of `test_input.sam` that can be used as input into Deduper
    - This file was created by running the `samtools sort` command from the top level directory of this Git repository:
        ```
        samtools sort -O SAM -o unit_tests/test_input_sorted.sam unit_tests/test_input.sam
        ```
- `expected_test_output.sam` contains Deduper's expected output SAM file after after running Deduper from the top level directory of this Git repository:
    ```
    ./karumuru_deduper.py -f unit_tests/test_input_sorted.sam -o unit_tests/expected_test_output.sam -u example_files/STL96.txt -s unit_tests/expected_test_output_stats.txt
    ```
- `expected_test_output_stats.txt` contains corresponding deduplication statistics calculated by Deduper


## Test Cases
The following table outlines the reads (in order) in `test_input_sorted.sam`, and explains why each read is retained or removed by Deduper: <br>
| Read Number | Chromosome | Strand | UMI | 5' Start Position | Comments
| - | - | - | - | - | - |
| 1 | 2 | + | AGAGTCCG | 52159545 | Not known UMI -> **Remove** |
| 2 | 2 | + | GGATAACG | 52308305 | Known UMI & 5' start position not encountered before -> Unique read -> **Retain** |
| 3 | 2 | + | CTAGGAAG | 52308307 - 2 = **52308305** *(accounting for soft clipping)* | Known UMI & 5' start position not encountered before -> Unique read -> **Retain** |
| 4 | 2 | + | GGATAACG | 52328864 | Known UMI & 5' start position not encountered before -> Unique read -> **Retain** |
| 5 | 2 | + | AGCTACCA | 76710746 | Known UMI & 5' start position not encountered before -> Unique read -> **Retain** |
| 6 | 2 | + | AGCTACCA | 76710746 | Same chromosome, strand, UMI, 5' start position as previous read -> Duplicate read -> **Remove** |
| 7 | 2 | + | TAGCAAGG | 76718924 | Known UMI & 5' start position not encountered before -> Unique read -> **Retain** |
| 8 | 2 | + | CTGTTCAC | 76814284 | Known UMI & 5' start position not encountered before -> Unique read -> **Retain** |
| 9 | 2 | + | AGACACTC | 76906303 | Known UMI & 5' start position not encountered before -> Unique read -> **Retain** |
| 10 | 2 | + | AGACACTC | 76906303 - 5 = **76906298** *(accounting for soft clipping)* | Known UMI & 5' start position not encountered before -> Unique read -> **Retain** |
| 11 | 2 | - | AGACACTC | 76906303 + 71 = **76906374** *(finding start position of read on minus strand)* | Known UMI & first read on minus strand encountered -> Unique read -> **Retain** |
| 12 | 2 | + | AGACACTC | 76906305 - 2 = **76906303** *(accounting for soft clipping)* | Same chromosome, strand, UMI, 5' start position as Read #9 -> Duplicate read -> **Remove** |
| 13 | 2 | - | AGACACTC | 76906306 + 46 + 2 + 20 = **76906374** *(finding start position of read on minus strand)* | Same strand, UMI, 5' strart position as Read #11 -> Duplicate read -> **Remove** |
| 14 | 3 | + | AGTGCTGT | 76776807 | Known UMI & first read encountered on this chromosome -> Unique read -> **Retain** |
| 15 | 3 | + | AGTGCTGT | 76776808 - 1 = **76776807** *(accounting for soft clipping)* | Same chromosome, strand, UMI, 5' start position as previous read -> Duplicate read -> **Remove** |
| 16 | 16 | + | AGCTACCA | 76710746 | Known UMI & first read encountered on this chromosome -> Unique read -> **Retain** |