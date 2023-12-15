import pathlib
import itertools

configfile:
    "config.yaml"

TRANSCRIPTOME_FASTA = config['reference_genome']['fasta']
GTF = config['reference_genome']['gtf']

ITERS_PER_BATCH = config['optimization']['iters_per_batch']
TAKE_TOP_N = config['optimization']['take_top_n']
NUM_BATCHES = config['optimization']['num_batches']
R1_FASTQ = config['real_sample']['FASTQ_R1']
R2_FASTQ = config['real_sample']['FASTQ_R2']
BAM = config['real_sample']['BAM']
BAM_PATH = pathlib.Path(BAM)
BAM_SORTED = BAM_PATH.parent / f"{BAM_PATH.stem}.sorted.bam"
BAI = BAM_PATH.parent / f"{BAM_PATH.stem}.sorted.bam.bai"

wildcard_constraints:
    sample = "[0-9]+", # actually only use sample=1 in this
    run = "[a-zA-Z0-9_]+", # these will actually be of the form batchX_XX in this, where X's represent digits
    filename = "[a-zA-Z0-0_\\.]+"

rule all:
    input:
        "results/compare_real_sim_cov",
        expand("data/batch{num_batches}_{num}/scores.json", num_batches=[NUM_BATCHES], num=range(ITERS_PER_BATCH)),

rule prep_input:
    input:
        gene_ids = "chosen_transcripts.txt",
        transcriptome = TRANSCRIPTOME_FASTA,
        gtf = GTF,
    output:
        molecules = "input_data/sample1/molecule_file.txt",
        reference = "input_data/reference_genome.fasta"
    params:
        num_molecules = 250
    resources:
        mem_mb = 12_000
    script:
        "scripts/prep_input.py"

rule generate_config:
    input:
        scores = lambda wildcards: [f"data/batch{batch}_{num}/scores.json"
                                        for batch in range(int(wildcards.batch))
                                        for num in range(ITERS_PER_BATCH if batch > 0 else 1)],
        configs = lambda wildcards: [f"config/generated/batch{batch}_{num}.json"
                                        for batch in range(int(wildcards.batch))
                                        for num in range(ITERS_PER_BATCH if batch > 0 else 1)],
    output:
        config = expand("config/generated/batch{{batch}}_{num}.json", num = range(ITERS_PER_BATCH))
    params:
        top_n = TAKE_TOP_N
    script:
        "scripts/generate_config.py"

rule run_beers:
    input:
        config_template = "config/config.template.yaml",
        molecule_files = "input_data/sample1/molecule_file.txt",
        reference_genome = "input_data/reference_genome.fasta",
        config = "config/generated/batch{batch}_{num}.json",
    output:
        flag_file = "data/batch{batch}_{num}/beers/finished_flag"
    params:
        molecule_dir = "input_data/",
        beers_dir = directory("data/batch{batch}_{num}/beers/"),
        config = "data/batch{batch}_{num}/beers.config.yaml",
    resources:
        mem_mb = 18_000
    run:
        #Generate config template
        import numpy as np
        import string
        import json
        config_fill_values = json.load(open(input.config))
        config_fill_values['camparee_output'] = "none"# input.camparee_output
        config_fill_values['molecule_files'] = pathlib.Path(params.molecule_dir).resolve()
        config_fill_values['reference_genome'] = pathlib.Path(input.reference_genome).resolve()
        # Since the frequencies sum to 1, we fill in the remaining T frequency from the other three
        config_fill_values['fwd_T_freq'] = list(1 - (np.array(config_fill_values['fwd_A_freq']) + np.array(config_fill_values['fwd_C_freq']) + np.array(config_fill_values['fwd_G_freq'])))
        config_fill_values['rev_T_freq'] = list(1 - (np.array(config_fill_values['rev_A_freq']) + np.array(config_fill_values['rev_C_freq']) + np.array(config_fill_values['rev_G_freq'])))

        config_template = string.Template(open(input.config_template, "r").read())
        config = config_template.substitute(config_fill_values)
        with open(params.config, "w") as f:
            f.write(config)
        # Run beers
        beers_cmd =  f"beers --configfile {params.config} --directory {params.beers_dir} -c 1 --rerun-incomplete"
        print(beers_cmd)
        shell(beers_cmd)
        shell("touch {output.flag_file}")

rule compute_gc_content_from_beers:
    input:
        "data/{run}/beers/finished_flag"
    output:
        gc_content = "data/{run}/sample{sample}/gc_content.txt",
    params:
        fastq_files = lambda wildcards: pathlib.Path(f"data/{wildcards.run}/beers/results/").glob(f"S{wildcards.sample}_L1*.fastq")
    resources:
        mem_mb = 6_000,
    script:
        "scripts/compute_gc_content.py"

rule compute_gc_content_from_real:
    input:
        R1_FASTQ,
        R2_FASTQ,
    output:
        gc_content = "real_data/gc_content.txt"
    params:
        fastq_files = [R1_FASTQ, R2_FASTQ]
    resources:
        mem_mb = 6_000,
    script:
        "scripts/compute_gc_content.py"

rule generate_BAM:
    input:
        "data/{run}/beers/finished_flag",
    output:
        bam = "data/{run}/sample{sample}/BEERS_output.sorted.bam",
    params:
        data_folder = "data/{run}/beers/results/",
    run:
        import uuid
        tmpdir = pathlib.Path(resources.tmpdir) / str(uuid.uuid4())
        tmpdir.mkdir(exist_ok=True)
        files = pathlib.Path(params.data_folder).glob(f"S{wildcards.sample}_L*.sam")
        print(f"Temp dir: {tmpdir}")
        for sam in files:
            print(f"Sorting on {sam}")
            shell(f"samtools sort -o {tmpdir}/{sam.name}.sorted.bam -O bam {sam}")
        shell(f"samtools merge {output.bam} {tmpdir}/*.sorted.bam")

rule compute_coverage:
    input:
        bam = "data/{run}/sample{sample}/BEERS_output.sorted.bam"
    output:
        "data/{run}/sample{sample}/coverage.txt"
    shell:
        "samtools depth -a -s -d 0 {input} > {output}"

rule summarize_coverage:
    input:
        cov = "{path}/coverage.txt",
    output:
        summary = "{path}/coverage_summary.txt"
    script:
        "scripts/summarize_coverage.py"

rule compute_real_coverage:
    input:
        bam = BAM_SORTED,
        gene_ids = "chosen_transcripts.txt",
        gtf = GTF,
        bai = BAI,
    output:
        cov = "real_data/coverage.txt"
    script:
        "scripts/compute_real_coverage.py"

rule compare_real_sim_cov:
    input:
        sim_full_cov = "data/batch5_14/sample1/coverage.txt",
        sim_cov = "data/batch5_14/sample1/coverage_summary.txt",
        sim_gc = "data/batch5_14/sample1/gc_content.txt",
        sim_seq = "data/batch5_14/sample1/seq_frequencies.json",
        sim_frag = "data/batch5_14/sample1/frag_sizes.txt",
        real_full_cov = "real_data/coverage.txt",
        real_cov = "real_data/coverage_summary.txt",
        real_gc = "real_data/gc_content.txt",
        real_seq = "real_data/seq_frequencies.json",
        real_frag = "real_data/frag_sizes.txt",
        reference_genome = "input_data/reference_genome.fasta",
    output:
        outdir = directory("results/compare_real_sim_cov/")
    script:
        "scripts/compare_real_sim_cov.py"

rule compute_seq_bias:
    input:
        flag_file = "data/{run}/beers/finished_flag",
    params:
        fastq_files = lambda wildcards: pathlib.Path(f"data/{wildcards.run}/beers/results/").glob(f"S{wildcards.sample}_L1*.fastq")
    output:
        seq_frequencies = "data/{run}/sample{sample}/seq_frequencies.json",
    script:
        "scripts/compute_seq_bias.py"

rule compute_seq_bias_from_real:
    input:
        R1_FASTQ,
        R2_FASTQ,
    params:
        fastq_files = [R1_FASTQ, R2_FASTQ]
    output:
        seq_frequencies = "real_data/seq_frequencies.json",
    script:
        "scripts/compute_seq_bias.py"

rule compute_frag_sizes_from_beers:
    input:
        bam = "data/{run}/sample1/BEERS_output.sorted.bam",
        gene_ids = "chosen_transcripts.txt",
        bai = "data/{run}/sample1/BEERS_output.sorted.bam.bai",
    output:
        frag_sizes = "data/{run}/sample1/frag_sizes.txt"
    script:
        "scripts/compute_frag_sizes.py"

rule compute_frag_sizes_from_real:
    input:
        bam = BAM_SORTED,
        gene_ids = "chosen_transcripts.txt",
        gtf = GTF,
        bai = BAI,
    output:
        frag_sizes = "real_data/frag_sizes.txt"
    script:
        "scripts/compute_frag_sizes_from_real.py"

rule compute_scores:
    input:
        sim_cov = "data/{run}/sample1/coverage_summary.txt",
        sim_gc = "data/{run}/sample1/gc_content.txt",
        sim_seq = "data/{run}/sample1/seq_frequencies.json",
        sim_frag = "data/{run}/sample1/frag_sizes.txt",
        real_cov = "real_data/coverage_summary.txt",
        real_gc = "real_data/gc_content.txt",
        real_seq = "real_data/seq_frequencies.json",
        real_frag = "real_data/frag_sizes.txt",
    output:
        results = "data/{run}/scores.json"
    script:
        "scripts/compute_scores.py"

rule sort_bam:
    input:
        bam = "{path}/{filename}.bam"
    output:
        sortedbam = "{path}/{filename}.sorted.bam"
    shell:
        "samtools sort -o {output.sortedbam} {input.bam} "

rule generate_BAI:
    input:
        bam = "{path}/{filename}.sorted.bam",
    output:
        bai = "{path}/{filename}.sorted.bam.bai",
    resources:
        mem_mb = 6_000
    shell:
        'samtools index -b {input} >> {output}'
