import polars as pl
import gzip
import gtfparse

NUM_MOLECULES = snakemake.params.num_molecules

with open(snakemake.input.gene_ids) as f:
    gene_ids = [l.strip() for l in f.read().splitlines()]

# Grab their transcript sequences
with gzip.open(snakemake.input.transcriptome, "rt") as f:
    transcriptome = {}
    transcript_id = None
    sequence = []
    for line in f:
        if line.startswith(">"):
            transcriptome[transcript_id] = ''.join(sequence)
            transcript_id = line.removeprefix(">").split(".")[0]
            sequence = []
        else:
            sequence.append(line.strip())
    transcriptome[transcript_id] = ''.join(sequence)


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

with open(snakemake.output.molecules, "wt") as molecules:
    for transcript in transcript_ids['transcript_id']:
        seq = transcriptome[transcript] + "A"*200 # polyA tail
        matchlen = len(seq)
        #transcript_id chrom parental_start parental_cigar ref_start ref_cigar strand sequence
        for i in range(NUM_MOLECULES):
            molecules.write(f"{transcript}\t{transcript}\t1\t{matchlen}M\t1\t{matchlen}M\t+\t{seq}\n")

with open(snakemake.output.reference, "wt") as reference:
    for transcript in transcript_ids['transcript_id']:
        seq = transcriptome[transcript]
        reference.write(f">{transcript}\n")
        reference.write(f"{seq}\n")
