import json
import numpy as np
import copy
import functools

## Everytime batch we reduce by half the perturbation we make
SCALE_FACTOR = 2**(1-int(snakemake.wildcards.batch))

# Load performance scores of all the previous rounds
scores = {}
for num, score in enumerate(snakemake.input.scores):
    with open(score) as scorefile:
        scores[num] = json.load(scorefile)['overall']

# Take the top N of them
# Lower scores are better
sorted_scores = sorted(scores.keys(), key = lambda num: scores[num])
best_scorers = sorted_scores[:snakemake.params.top_n]
print(best_scorers)

@functools.cache
def get_config(num):
    with open(snakemake.input.configs[num]) as base_config_file:
        base_config = json.load(base_config_file)
    return base_config

# Generate a random new iteration from the previous one
# first load one of the top scorers from the previous round to initialize this
rng = np.random.default_rng([snakemake.wildcards.batch])
for outfile in snakemake.output.config:
    ## Choose the base configuration which we modify
    initialize_from = rng.choice(best_scorers)
    base_config = get_config(initialize_from)

    # Perturb the base config to get the new config
    config = copy.copy(base_config)
    def perturb_num(name, std, mode="linear", min_val=None, max_val=None):
        if mode == "linear":
            config[name] = rng.normal(config[name], std*SCALE_FACTOR)
        elif mode == "log":
            config[name] = np.exp(rng.normal(np.log(config[name]), std*SCALE_FACTOR))
        if min_val is not None:
            config[name] = max(config[name], min_val)
        if max_val is not None:
            config[name] = min(config[name], max_val)
    perturb_num("gc_bias_constant", 0.1)
    perturb_num("gc_bias_linear", 0.1)
    perturb_num("gc_bias_quadratic", 10)
    perturb_num("breakpoint_prob_per_base", 1, mode="log", max_val=1)
    perturb_num("primes_per_kb", 5, min_val=1)
    perturb_num("fragmentation_rate", 1, mode="log", max_val=1)

    # Ensure that these four values are in increasing order by perturbing their differences
    config['mid_min_diff'] = config['frag_mid_start_size'] - config['frag_min_size']
    config['mid_mid_diff'] = config['frag_mid_end_size'] - config['frag_mid_start_size']
    config['max_mid_diff'] = config['frag_max_size'] - config['frag_mid_end_size']
    perturb_num("frag_min_size", 30, min_val = 50)
    perturb_num("mid_min_diff", 30, min_val=1)
    perturb_num("mid_mid_diff", 30, min_val=1)
    perturb_num("max_mid_diff", 30, min_val=1)
    config['frag_mid_start_size'] = config['frag_min_size'] + config['mid_min_diff']
    config['frag_mid_end_size'] = config['frag_mid_start_size'] + config['mid_mid_diff']
    config['frag_max_size'] = config['frag_mid_end_size'] + config['max_mid_diff']

    # NOTE: for now, we don't vary the sequence frequencies, and just use the empirical frequencies
    #for base in ["A", "C", "G"]:
    #    for dir in ['fwd', 'rev']:
    #        current = np.array(config[f"{dir}_{base}_freq"])


    with open(outfile, "wt") as out:
        json.dump(config, out)
