import polars as pl
import scipy.stats
import numpy as np
import json

def EMD(a,b):
    # Earth-movers distance between two distributions
    # note: use scipy.stats.wasserstein_distance to compute this
    #       when you just have two samples, not distributions
    assert len(a) == len(b)
    # Normalize to a probability distribution
    a = np.array(a, dtype=float)
    a /= a.sum()
    b = np.array(b, dtype=float)
    b /= b.sum()

    emd = np.zeros(len(a))
    for i in range(1,len(a)):
        emd[i] = a[i-1] + emd[i-1] - b[i-1]

    # and we normalize the end result to b's distribution for comparison
    # across variables
    # Normalize by the standard deviation (computed from the distribution b)
    b_val = np.arange(len(b))
    b_mean = (b_val*b).sum()
    b_sd = np.sqrt((((b_val - b_mean)**2)*b).sum())
    return np.abs(emd).sum() / b_sd

### COVERAGE
# scored by earth-movers distance normalized by SD in real data
sim_cov = pl.read_csv(snakemake.input.sim_cov, separator="\t")
real_cov = pl.read_csv(snakemake.input.real_cov, separator="\t")

cov_dist = scipy.stats.wasserstein_distance(
        sim_cov.filter(pl.col("transcript_length") < 1500)['coef_of_var'],
        real_cov.filter(pl.col("transcript_length") < 1500)['coef_of_var'],
    ) / real_cov.filter(pl.col("transcript_length") < 1500)['coef_of_var'].std()
depth_dist = scipy.stats.wasserstein_distance(
        sim_cov.filter(pl.col("transcript_length") > 1500)['exp_depth_pos_regression'],
        real_cov.filter(pl.col("transcript_length") > 1500)['exp_depth_pos_regression'],
    ) / real_cov.filter(pl.col("transcript_length") > 1500)['exp_depth_pos_regression'].std()

### GC CONTENT
# scored by earth-movers distance normalized by SD in real data
sim_gc = pl.read_csv(snakemake.input.sim_gc, separator="\t")
real_gc = pl.read_csv(snakemake.input.real_gc, separator="\t")

gc_dist = EMD(
        sim_gc['read_count'],
        real_gc['read_count'],
    )


### SEQUENCE BIAS
# Scored by average base-pair frequency mismatch
sim_seq = json.load(open(snakemake.input.sim_seq))
real_seq = json.load(open(snakemake.input.real_seq))
def seq_to_mat(seq):
    return np.array([
        seq['fwd_frequencies']['A'],
        seq['fwd_frequencies']['C'],
        seq['fwd_frequencies']['G'],
        seq['fwd_frequencies']['T'],
        seq['rev_frequencies']['A'],
        seq['rev_frequencies']['C'],
        seq['rev_frequencies']['G'],
        seq['rev_frequencies']['T'],
    ])

seq_dist = np.abs(seq_to_mat(sim_seq) - seq_to_mat(real_seq)).mean()


### FRAG SIZES
sim_frag = pl.read_csv(snakemake.input.sim_frag, separator="\t") \
        .sort("frag_size")
real_frag = pl.read_csv(snakemake.input.real_frag, separator="\t") \
        .sort("frag_size")
longest_frag = max(sim_frag['frag_size'].max(), real_frag['frag_size'].max())
all_frag_sizes = pl.DataFrame({"frag_size": np.arange(1, longest_frag+1)})
sim_frag = sim_frag.join(all_frag_sizes, on="frag_size", how="outer") \
        .with_columns(pl.col("counts").fill_null(0).alias("counts"))
real_frag = real_frag.join(all_frag_sizes, on="frag_size", how="outer") \
        .with_columns(pl.col("counts").fill_null(0).alias("counts"))

frag_dist = EMD(sim_frag['counts'], real_frag['counts'])

## COMBINE THEM FOR AN OVERALL SCORE
overall = cov_dist + depth_dist + gc_dist + seq_dist + frag_dist

## OUTPUT RESULTS
results = {
    "cov_dist": cov_dist,
    "depth_dist": depth_dist,
    "gc_dist": gc_dist,
    "seq_dist": seq_dist,
    "frag_dist": frag_dist,
    "overall": overall
}
with open(snakemake.output.results, "wt") as f:
    json.dump(results, f)
