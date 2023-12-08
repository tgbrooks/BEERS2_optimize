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
real = pl.read_csv("real_data/WT4_PolyA/coverage_summary.txt", separator="\t") \
        .filter(pl.col("mean") > 100) \
        .with_columns(pl.lit("real").alias("source"))

sim = pl.read_csv("data/all_bias/sample1/coverage_summary.txt", separator="\t") \
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
