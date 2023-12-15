#!/usr/bin/env sh
set -e
module load /project/itmatlab/sharedmodules/use.shared
module load samtools/1.11
module load htslib/1.11
module load /project/itmatlab/sharedmodules/OpenSSL-v3.0.1
source ~/BEERS2_assess/venv311/bin/activate


### SETUP NOTE:
# In order to use the following --profile lsf command
# you need to follow the instructions at https://github.com/Snakemake-Profiles/lsf
# and set up LSF support for snakemake

bsub -e logs/snakemake.err \
     -o logs/snakemake.out \
     snakemake --profile lsf -j 50 -c 50 "$@"
