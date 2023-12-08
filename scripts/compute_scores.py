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
    a = np.array(a)
    a /= a.sum()
    b = np.array(b)
    b /= b.sum()

    emd = np.zeros(len(a))
    for i in range(1,len(a)):
        emd[i] = a[i-1] + emd[i-1] - b[i-1]
    return np.abs(emd).sum()


### COVERAGE
sim_cov = pl.read_csv(snakemake.input.sim_cov, separator="\t")
real_cov = pl.read_csv(snakemake.input.real_cov, separator="\t")

cov_dist = scipy.stats.wasserstein_distance(
        sim_cov['coef_of_var'],
        real_cov['coef_of_var'],
    )
depth_dist = scipy.stats.wasserstein_distance(
        sim_cov['depth_pos_regression'] / sim_cov['mean'],
        real_cov['depth_pos_regression'] / real_cov['mean'],
    )

### GC CONTENT
sim_gc = pl.read_csv(snakemake.input.sim_gc, separator="\t")
real_gc = pl.read_csv(snakemake.input.real_gc, separator="\t")

gc_dist = EMD(
        sim_gc['read_count'],
        real_gc['read_count'],
    )


### SEQUENCE BIAS
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

seq_dist = np.abs(seq_to_mat(sim_seq) - seq_to_mat(real_seq)).sum()

results = {
    "cov_dist": cov_dist,
    "depth_dist": depth_dist,
    "gc_dist": gc_dist,
    "seq_dist": seq_dist,
}
with open(snakemake.output.results, "wt") as f:
    json.dump(results, f)
