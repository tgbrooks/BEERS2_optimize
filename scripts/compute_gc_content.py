import pathlib
import re

import pandas
import numpy as np
import gzip

N_SAMPLES = 100_000 # Per fastq

def maybe_gzip_open(path):
    if str(path).endswith("gz"):
        return gzip.open(path, "rt")
    return open(path)

bin_edges = np.linspace(0,1,101)
counts_by_bin = np.zeros(100)
for fastq_file in snakemake.params.fastq_files:
    print("Processing", fastq_file)
    with maybe_gzip_open(fastq_file) as fastq:
        id = None
        i = 0
        while True:
            line = fastq.readline()
            if line == '':
                break
            i += 1
            assert line.startswith("@"), f"Error in {fastq_file} on line {i}: expected @ not found instead received {repr(line)}"
            seq = fastq.readline().strip()
            GC_content = (seq.count('C') + seq.count('G') + seq.count('N')/2) / len(seq)
            bin = int(GC_content*100)
            counts_by_bin[min(bin, 99)] += 1
            fastq.readline() # + -- skip
            fastq.readline() # quality -- skip
            if i > N_SAMPLES:
                break

GC_content = pandas.DataFrame({
        "lower_gc_content": bin_edges[:-1],
        "upper_gc_content": bin_edges[1:],
        "read_count": counts_by_bin,
    })
GC_content.to_csv(snakemake.output.gc_content, sep="\t", index=False)
#GC_by_transcript = GC_content.groupby(["run", "sample", "TranscriptID"]).GC_content.mean()
#GC_by_transcript.reset_index().to_csv(snakemake.output.gc_content, sep="\t")
