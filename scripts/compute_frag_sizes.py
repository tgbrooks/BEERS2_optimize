import polars as pl
import pandas
import numpy as np
import io
import gzip
import gtfparse
import subprocess
import beers_utils.cigar

## Computes frag sizes of our simulated data
N_SAMPLES_PER_TRANSCRIPT = 1_000

bam_path = snakemake.input.bam

with open(snakemake.input.gene_ids) as f:
    gene_ids = [l.strip() for l in f.read().splitlines()]

frag_sizes = []
def get_match_len(cigar):
    return beers_utils.cigar.match_seq_length(beers_utils.cigar.split_cigar(cigar))

# Read the bam file from this location
cmd = ["samtools", "view", bam_path]
print(' '.join(cmd))
p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

working_set = dict()
for line in iter(p.stdout.readline, b''):
    if line == '':
        break
    line = line.decode("utf-8")
    id, flag_, _, pos_, _, cigar, *_other = line.split('\t')
    flag = int(flag_)
    pos = int(pos_)

    if flag & 0x100:
        # Secondary alignment - don't want
        continue

    if id in working_set:
        (mate_pos, mate_cigar) = working_set[id]
        if mate_pos < pos:
            start = mate_pos
            end = pos + get_match_len(cigar)
        else:
            start = pos
            end = mate_pos + get_match_len(mate_cigar)

        frag_sizes.append(end - start + 1)
        del working_set[id]
    else:
        working_set[id] = (pos, cigar)

all_frag_sizes = pl.DataFrame({"frag_size": frag_sizes})['frag_size'].value_counts()

all_frag_sizes.write_csv(snakemake.output.frag_sizes, separator="\t", has_header=True)
