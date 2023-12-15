import pathlib
import polars as pl
import pandas
import seaborn as sns
import numpy as np

sim_dir = pathlib.Path(snakemake.input.sim)

outdir = pathlib.Path(snakemake.output.outdir)
outdir.mkdir(exist_ok=True)

order = ["real", "sim"]

# Load gene lengths
with open(snakemake.input.reference_genome, "rt") as f:
    transcript_lengths = {}
    while True:
        line = f.readline().strip()
        if line == '':
            break
        transcript_id = line.removeprefix(">").strip()
        seq = f.readline().strip()
        transcript_lengths[transcript_id] = len(seq)

real = pl.read_csv(snakemake.input.real_cov, separator="\t") \
        .with_columns(pl.lit("real").alias("source"))

sim = pl.read_csv(sim_dir / "coverage_summary.txt", separator="\t") \
        .with_columns(pl.lit("sim").alias("source"))

both = pl.concat([
        real,
        sim
    ]) \
    .with_columns( (pl.col("depth_pos_regression") / pl.col("mean")).alias("regression") ) \
    .with_columns(pl.col("transcript_id").map_dict(transcript_lengths).alias("transcript_length")) \
    .to_pandas()
both['length'] = pandas.cut(both['transcript_length'], [0,3000, max(transcript_lengths.values())+1])
print(both)

fig = sns.displot(
    data = both,
    hue = "source",
    x = "coef_of_var",
    hue_order = order,
    kind = "kde",
)
fig.savefig(outdir / "coef_of_var.png", dpi=300)

fig = sns.displot(
    both,
    x = "coef_of_var",
    hue = "source",
    hue_order = order,
    col =  "length",
    kind = "kde",
    facet_kws=dict(sharex=False, sharey=False),
)
fig.savefig(outdir / "coef_of_var.by_transcript_length.png", dpi=300)

fig = sns.displot(
    data = both,
    hue = "source",
    x = "depth_pos_corr",
    hue_order = order,
    kind = "kde",
)
fig.savefig(outdir / "depth_pos_corr.png", dpi=300)

fig = sns.displot(
    data = both,
    hue = "source",
    x = "regression",
    hue_order = order,
    kind = "kde",
)
fig.savefig(outdir / "regression.png", dpi=300)

fig = sns.displot(
    both,
    x = "regression",
    hue = "source",
    hue_order = order,
    col =  "length",
    kind = "kde",
)
fig.savefig(outdir / "regression.by_transcript_length.png", dpi=300)

fig = sns.displot(
    both,
    x = "exp_depth_pos_regression",
    hue = "source",
    hue_order = order,
    col =  "length",
    kind = "kde",
    facet_kws=dict(sharex=False, sharey=False),
)
fig.savefig(outdir / "exp_regression.by_transcript_length.png", dpi=300)



### GC CONTENT
sim_gc = pl.read_csv(sim_dir / "gc_content.txt", separator="\t") \
        .with_columns(
                (pl.col('read_count') / pl.col('read_count').mean()).alias('density'),
                # Summarize to the nearest 5% - the plot is ragged otherwise due to floating point issues
                ((pl.col('lower_gc_content') * 20).cast(int) / 20 + 0.025).alias('GC'),
                pl.lit("sim").alias("source")
            ) \
        .groupby("GC") \
        .agg(pl.col('density').sum(), pl.col("source").first())
real_gc = pl.read_csv(snakemake.input.real_gc, separator="\t") \
        .with_columns(
                (pl.col('read_count') / pl.col('read_count').mean()).alias('density'),
                # Summarize to the nearest 5% - the plot is ragged otherwise due to floating point issues
                ((pl.col('lower_gc_content') * 20).cast(int) / 20 + 0.025).alias('GC'),
                pl.lit("real").alias("source")
            ) \
        .groupby("GC") \
        .agg(pl.col('density').sum(), pl.col("source").first())
both_gc = pl.concat([
        real_gc,
        sim_gc
    ]).to_pandas()

fig = sns.relplot(
    both_gc,
    x = "GC",
    y = "density",
    hue = "source",
    hue_order = order,
    kind = "line",
)
fig.savefig(outdir / "gc.png", dpi=300)

### FRAG SIZES
sim_frag = pl.read_csv(sim_dir / "frag_sizes.txt", separator="\t") \
        .sort("frag_size") \
        .with_columns(pl.lit("sim").alias("source"))
real_frag = pl.read_csv(snakemake.input.real_frag, separator="\t") \
        .sort("frag_size")\
        .with_columns(pl.lit("real").alias("source"))

both_frag = pl.concat([
        sim_frag,
        real_frag,
    ]) \
    .with_columns((pl.col('counts') / pl.col('counts').sum().over('source')).alias("density"))\
    .to_pandas()

fig = sns.relplot(
    both_frag,
    x = "frag_size",
    y = "density",
    hue = "source",
    hue_order = order,
    kind = "line",
)
fig.savefig(outdir / "frag_size.png", dpi=300)
