from read import read_fast5
from normalization import med_mad, normalize_signal
import numpy as np
from matplotlib import pyplot as plt

read_data = read_fast5('data/0a0bdc5c-8f8f-41ea-a4d1-4ff6344fac3e.fast5')
read_data = read_data[list(read_data.keys())[0]]

# check if it has a segmentation table
if read_data.segmentation is None:
    print('This file does not have a reference')

# normalize the raw signal
med, mad = med_mad(read_data.raw, factor=1.0)
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

base_pos = np.where(dna_in_raw != '')[0]
dna_in_raw_base_idx = np.full(dna_in_raw.shape, -1)
dna_in_raw_base_idx[base_pos] = np.arange(len(base_pos))

segment_st = 1000
segment_nd = 1400

plt.figure(figsize=(10, 5))
plt.plot(norm_signal[segment_st:segment_nd], color = 'black')

for b in dna_in_raw_base_idx[segment_st:segment_nd]:
    if b == -1:
        continue

    st = dna_bases_start_in_raw[b] - segment_st
    nd = dna_bases_end_in_raw[b] - segment_st
    
    plt.plot(
        [st, nd],
        [read_data.segmentation['norm_mean'][b], read_data.segmentation['norm_mean'][b]],
        color = 'red'
    )

    plt.text(
        x = st + ((nd-st)/2),
        y = read_data.segmentation['norm_mean'][b] + 0.2,
        s = read_data.segmentation['base'][b],
        color = 'red'
    )