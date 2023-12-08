import pathlib
import numpy as np
import gzip
import json

N_SAMPLES = 100_000 # Per fastq
BASES = ['A', 'C', 'G', 'T']
PRIMER_LENGTH = 6 
N_SAMPLES = 100_000

def maybe_gzip_open(path):
    if str(path).endswith("gz"):
        return gzip.open(path, "rt")
    return open(path)

fwd_count_matrix = np.zeros((4, PRIMER_LENGTH))
rev_count_matrix = np.zeros((4, PRIMER_LENGTH))
reverse = False # We assume that the reverse reads are the second fastq file
for fastq_file in snakemake.params.fastq_files:
    count_matrix = fwd_count_matrix if not reverse else rev_count_matrix
    print("Processing", fastq_file)
    with maybe_gzip_open(fastq_file) as fastq:
        i = 0
        while True:
            line = fastq.readline()
            if line == '':
                break
            i += 1
            assert line.startswith("@"), f"Error in {fastq_file} on line {i}: expected @ not found instead received {repr(line)}"

            seq = fastq.readline().strip()
            count_matrix += [[1 if x == char else 0 for x in seq[:PRIMER_LENGTH]] for char in BASES]

            fastq.readline() # + -- skip
            fastq.readline() # quality -- skip
            if i > N_SAMPLES:
                break
    reverse = True

num_molecules = i
results = { 
        # Frequencies of bases by position for the forward read
        "fwd_frequencies": {base: list(fwd_count_matrix[i] / num_molecules) for i, base in enumerate(BASES)},
        # Frequencies of bases by position for the reverse read
        "rev_frequencies": {base: list(rev_count_matrix[i] / num_molecules) for i, base in enumerate(BASES)},
}
json.dump(results, open(snakemake.output.seq_frequencies, 'w'))
