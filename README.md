# BEERS2 Parameter Optimization

When using BEERS2 to mimic a specific real sample, this pipeline can help you to determine which parameter settings to use.
This works by a simple 'particle swarm' optimization approach where batches of 50 BEERS2 runs on small input samples are performed and evaluated compared to a real sample.
The best performing parameters are then selected to seed the next batch, with each iteration using a perturbation of the seed parameters.
Each batch reduces the size of the perturbation compared to the previous batch so as to be able to narrow in to high-performing parameters.

## SETUP

First, make sure that you have `samtools` installed and on the PATH.
We used `samtools` version 1.11.
Similarly, we require `python` version 3.11.

``` sh
git clone gitt@github.com:tgbrooks/BEERS2_optimization
cd BEERS2_optimization
python -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Configuring the optimization

Edit `config.yaml` to point to your reference genome files and the locations of the FASTQ and BAM files for the sample you want to mimic.

Then edit the `chosen_transcripts.txt` to contain a list of gene (not transcript) identifiers that you wish to use as a reference for the optimization process.
We recommend choosing genes that have high expression in the reference sample and have only a single annotated isoform (when using the 'basic' annotation from ENSEMBL).
The identifiers must match those in the GTF file.

By default, GC bias, fragmentation, size selection, and hexamer priming for strand synthesis are optimized.
If you wish to optimize other BEERS2 parameters, that can be done by modifying `config/config.template.yaml` to replace a parameter value by `$my_variable`, then editing `scripts/generate_config.py` to control how values of this are perturbed between batches.
You will also need to include a value for `my_variable` in the initial configuration, config/generated/batch0_0.json`, and you may need to modify the `run_beers` step in `Snakefile` if any special handling of the value is required.

## Running the optimization

The optimization can be started with a specific number of cores (for example 10) by running:

``` sh
snakemake -c 10
```

Alternateively, we recommend using a `snakemake` profile so that `snakemake` can make use of a cluster.

## Results

The optimized BEERS2 configuration is in `results/optimal_configuration.yaml`.
Figures comparing the final optimized configuration to the real sample are in `results/compare_real_sim_cov`.
