import pathlib
import polars as pl
import pandas
import seaborn as sns

outdir = pathlib.Path(snakemake.output.outdir)
outdir.mkdir(exist_ok=True)

## Note we filter out a couple transcripts that had very little coverage
## We selected these transcripts from the PORT-derived quantifications which must have differed from
## the analysis here (different ENSEMBL version) and so the 'high expressed' genes from that quantification aren't
## high expressed here
real = pl.read_csv(snakemake.input.real_cov, separator="\t") \
        .filter(pl.col("mean") > 100) \
        .with_columns(pl.lit("real").alias("source"))

sim = pl.read_csv(snakemake.input.sim_cov, separator="\t") \
        .with_columns(pl.lit("sim").alias("source"))


both = pl.concat([
        real,
        sim
    ]) \
    .with_columns(
        (pl.col("depth_pos_regression") / pl.col("mean")).alias("regression")
    ).to_pandas()

fig = sns.catplot(
    data = both,
    x = "source",
    y = "coef_of_var",
    kind = "violin",
)
fig.savefig(outdir / "coef_of_var.png", dpi=300)

fig = sns.catplot(
    data = both,
    x = "source",
    y = "depth_pos_corr",
    kind = "violin",
)
fig.savefig(outdir / "depth_pos_corr.png", dpi=300)

fig = sns.catplot(
    data = both,
    x = "source",
    y = "regression",
    kind = "violin",
)
fig.savefig(outdir / "regerssion.png", dpi=300)



### GC CONTENT
sim_gc = pl.read_csv(snakemake.input.sim_gc, separator="\t") \
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
    kind = "line"
)
fig.savefig(outdir / "gc.png", dpi=300)
