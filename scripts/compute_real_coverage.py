import polars as pl
import pandas
import numpy as np
import io
import gzip
import gtfparse
import subprocess


bam_path = snakemake.input.bam

with open(snakemake.input.gene_ids) as f:
    gene_ids = [l.strip() for l in f.read().splitlines()]

gtf = gtfparse.read_gtf(snakemake.input.gtf)

# We use just the BASIC transcript annotations
tx_to_gene = gtf.filter(pl.col('transcript_id') != '') \
            .filter(pl.col("tag").str.contains("basic")) \
            .groupby("transcript_id") \
            .agg(pl.col('gene_id').first().alias("gene_id"))

# Every gene maps to exactly one transcritp
transcript_ids = pl.DataFrame({"gene_id": gene_ids}).join(
        tx_to_gene,
        on="gene_id",
        how="left"
    )

sequences = transcript_ids.join(gtf.filter(pl.col("feature") == "exon"), on="transcript_id", how="left") \
        .sort(["transcript_id", "start"]) \
        .select("transcript_id", "seqname", "start", "end", "strand")

all_covs = []
for (transcript_id, exons) in sequences.groupby("transcript_id"):
    pos = 1
    for (_, chrom, start, end, strand) in exons.iter_rows():
        # compute coverage of this exon
        exon_len = end - start + 1
        region_str = f"{chrom}:{start}-{end}"
        cmd = ["samtools", "depth", bam_path, "-r", region_str, '-s', '-d', '0', '-a']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        results = p.stdout.read().decode('utf-8')
        if results == '':
            continue
        raw_cov = pandas.read_csv(io.StringIO(results), sep="\\s+", header=None)
        raw_cov.columns = ['chr', 'pos', 'depth']
        direction = 1 if strand == "+" else -1

        # convert to transcript-based digits
        cov = pl.DataFrame({
            "transcript_id": transcript_id,
            "pos": pos + np.arange(exon_len),
            "depth": raw_cov['depth'].values[::direction],
        })
        all_covs.append(cov)
        pos += exon_len
all_covs = pl.concat(all_covs)

all_covs.write_csv(snakemake.output.cov, separator="\t", has_header=False)
