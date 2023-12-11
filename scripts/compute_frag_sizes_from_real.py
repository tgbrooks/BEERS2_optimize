import polars as pl
import pandas
import numpy as np
import io
import gzip
import gtfparse
import subprocess
import beers_utils.cigar

N_SAMPLES_PER_TRANSCRIPT = 1_000

bam_path = snakemake.input.bam

with open(snakemake.input.gene_ids) as f:
    gene_ids = [l.strip() for l in f.read().splitlines()]

gtf = gtfparse.read_gtf(snakemake.input.gtf)

# We use just the BASIC transcript annotations
tx_to_gene = gtf.filter(pl.col('transcript_id') != '') \
            .filter(pl.col("tag").str.contains("basic")) \
            .groupby("transcript_id") \
            .agg(pl.col('gene_id').first().alias("gene_id"))

# Every gene maps to exactly one transcript
transcript_ids = pl.DataFrame({"gene_id": gene_ids}).join(
        tx_to_gene,
        on="gene_id",
        how="left"
    )

transcripts = transcript_ids.join(gtf.filter(pl.col("feature") == "transcript"), on="transcript_id", how="left") \
        .filter(pl.col("tag").str.contains("basic")) \
        .sort(["transcript_id"]) \
        .select("transcript_id", "seqname", "start", "end")
exons = transcript_ids.join(gtf.filter(pl.col("feature") == "exon"), on="transcript_id", how="left") \
        .sort(["transcript_id", "start"]) \
        .select("transcript_id", "seqname", "start", "end")
        #.join(transcripts, on=["transcript_id", "seqname"], how="left")

frag_sizes = []
for (transcript_id, chrom, trans_start, trans_end) in transcripts.iter_rows():

    starts = []
    ends = []
    for (_, chrom, start, end) in exons.filter(pl.col("transcript_id") == transcript_id).iter_rows():
        starts.append(start)
        ends.append(end)
    starts = np.array(starts)
    ends = np.array(ends)
    exon_lengths = ends - starts + 1 # inclusive regions in GTFs
    cumlen = np.concatenate([[0], np.cumsum(exon_lengths)])
    def get_relative_position(pos):
        i = np.searchsorted(starts, pos, side="right") - 1
        return pos - starts[i] + cumlen[i]
    def get_match_len(cigar):
        return beers_utils.cigar.match_seq_length(cigar)
    def is_compatible(pos, cigar):
        # computes whether an alignment could be from this transcript
        curr_pos = pos
        for op, length in cigar:
            if beers_utils.cigar.consumes[op]['ref']:
                if ((starts <= curr_pos) & (ends >= curr_pos + length)).any():
                    curr_pos += length
                else:
                    return False
        return True

    # Read the bam file from this location
    region_str = f"{chrom}:{trans_start}-{trans_end}"
    cmd = ["samtools", "view", bam_path, region_str]
    print(' '.join(cmd))
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)

    num_samples = 0
    working_set = dict()
    for line in iter(p.stdout.readline, b''):
        line = line.decode("utf-8")
        id, flag_, _, pos_, _, cigar, *_other = line.split('\t')
        flag = int(flag_)
        pos = int(pos_)
        cigar = beers_utils.cigar.split_cigar(cigar)

        if flag & 0x100:
            # Secondary alignment - don't want
            continue
        if not is_compatible(pos, cigar):
            continue

        if id in working_set:
            (mate_pos, mate_cigar) = working_set[id]
            if mate_pos < pos:
                start = get_relative_position(mate_pos)
                end = get_relative_position(pos) + get_match_len(cigar)
            else:
                start = get_relative_position(pos)
                end = get_relative_position(mate_pos) + get_match_len(mate_cigar)
            frag_sizes.append(end - start + 1)
            num_samples += 1
            del working_set[id]
        else:
            working_set[id] = (pos, cigar)

        if num_samples > N_SAMPLES_PER_TRANSCRIPT:
            break

all_frag_sizes = pl.DataFrame({"frag_size": frag_sizes})['frag_size'].value_counts()

all_frag_sizes.write_csv(snakemake.output.frag_sizes, separator="\t", has_header=True)
