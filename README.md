# Nanopore data utilities

## Source

`read.py`: read functions for `.fast5`, `.fasta`, `.fastq` and `.fna` files.
`normalization.py`: several raw data normalization functions.

## Dependencies

- `numpy`
- `scipy`
- `ont_fast5_api`

## Demo

```
from read import read_fast5
from normalization import med_mad, normalize_signal
import numpy as np

read_data = read_fast5('data/0a0bdc5c-8f8f-41ea-a4d1-4ff6344fac3e.fast5')
read_data = read_data[list(read_data.keys())[0]]

# check if it has a segmentation table
if not read_data.segmentation:
    print('This file does not have a reference')

# normalize the raw signal
med, mad = med_mad(read_data.raw)
norm_signal = normalize_signal(read_data.raw, med, mad)

# read_data.segmentation contains the mapping between the DNA sequence and 
# the raw signal, but this is relative to the start of the segmentation
dna_bases = read_data.segmentation['base']
dna_bases_start_in_raw = read_data.segmentation['start'] + read_data.start_rel_to_raw
dna_bases_end_in_raw = dna_bases_start_in_raw + read_data.segmentation['length']

# here we create an array, the length of the raw signal, and annotate where
# each DNA base corresponds to the DNA signal
dna_in_raw = np.full(norm_signal.shape, '', dtype='U1')
dna_in_raw[dna_bases_start_in_raw] = dna_bases

```