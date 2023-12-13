import polars as pl
import numpy as np

df = pl.read_csv(snakemake.input.cov, sep="\t", has_header=False, new_columns=["transcript_id", "pos", "depth"])

def regress(params, trim_3prime = 1000):
    a,b = params
    A = np.array([
        np.ones_like(a),
        a
    ]).T
    A = A[:-trim_3prime]
    b = b[:-trim_3prime]
    coef, resid, rank, s  = np.linalg.lstsq(A, b, rcond=None)
    return coef[1]

def exp_regress(params, trim_3prime = 1000):
    a,b = params
    A = np.array([
        np.ones_like(a),
        a
    ]).T
    A = A[:-trim_3prime]
    b = b[:-trim_3prime]
    coef, resid, rank, s  = np.linalg.lstsq(A, np.log(b+1), rcond=None)
    return coef[1]

coef_var = df.groupby("transcript_id") \
        .agg(
            (pl.col("depth").std()  / pl.col("depth").mean()).alias("coef_of_var"),
            pl.col("depth").std().alias("std"),
            pl.col("depth").mean().alias("mean"),
            pl.pearson_corr("depth", "pos").alias("depth_pos_corr"),
            pl.apply(exprs=["pos", "depth"], function=regress).alias("depth_pos_regression"),
            pl.apply(exprs=["pos", "depth"], function=exp_regress).alias("exp_depth_pos_regression"),
            pl.col("pos").max().alias("transcript_length"),
        )

coef_var.write_csv(snakemake.output.summary, sep="\t")
