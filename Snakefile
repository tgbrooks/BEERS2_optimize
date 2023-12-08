import pathlib
import itertools

configfile:
    "config.yaml"

run_configs = config['run_configs']
sample_ids = config['sample_ids']
lanes_used = config['lanes_used']

RCLONE_REMOTE = "aws_igv_data" # Must first run `rclone config` and set up a remote with this name for uploading trackhubs to
BUCKET_NAME = "itmat.igv.data" # Bucket to upload to with rclone for trackhubs
BUCKET_DIR = f"BEERS2_ASSESS/DATA"

RUN_BEERS_PATH = 'beers_venv/bin/run_beers'
CAMPAREE_CONFIG = 'config/camparee_config.yaml'

TRANSCRIPTOME_FASTA = "/project/itmatlab/index/SALMON-1.9.0_indexes/GRCm38.ensemblv93/Mus_musculus.GRCm38.Ensembl.v93.cdna_ncrna_dna.fa.gz"
GTF = "/project/itmatlab/index/SALMON-1.9.0_indexes/GRCm38.ensemblv93/Mus_musculus.GRCm38.93.gtf.gz"

wildcard_constraints:
    sample = "[0-9]+",
    run = "[a-zA-Z0-9_]+",
    filename = "[a-zA-Z0-0_\\.]+"

rule all:
    input:
        #"results/browser_tracks/url.txt",
        #"results/gc_content/",
        #"results/coverage",
        #"results/seq_bias/fwd_seq_frequencies.png",
        #"results/frag_sizes/frag_sizes.fragment_counts.A.png",
        #"results/polya_tail/polya_tail_length.png",
        #"results/pcr_dupes/pcr_dupes.png",
        #"results/igv/url.txt",
        expand("data/{run}/beers/finished_flag", run = run_configs.keys()),
        expand("data/{run}/sample1/BEERS_output.sorted.bam", run = run_configs.keys()),
        expand("data/{run}/sample1/coverage_summary.txt", run = run_configs.keys()),
        expand("data/{run}/sample1/gc_content.txt", run = run_configs.keys()),
        "results/compare_real_sim_cov",
        "real_data/WT4_PolyA/coverage_summary.txt",
        "real_data/WT4_PolyA/gc_content.txt",

rule prep_input:
    input:
        gene_ids = "chosen_transcripts.txt",
        transcriptome = TRANSCRIPTOME_FASTA,
        gtf = GTF,
    output:
        molecules = "input_data/sample1/molecule_file.txt",
        reference = "input_data/reference_genome.fasta"
    params:
        num_molecules = 150
    resources:
        mem_mb = 12_000
    script:
        "scripts/prep_input.py"

rule run_beers:
    input:
        config_template = "config/config.template.yaml",
        molecule_files = "input_data/sample1/molecule_file.txt",
        reference_genome = "input_data/reference_genome.fasta",
    output:
        flag_file = "data/{run}/beers/finished_flag"
    params:
        molecule_dir = "input_data/",
        beers_dir = directory("data/{run}/beers/"),
        config = "data/{run}/beers.config.yaml",
    resources:
        mem_mb = 18_000
    run:
        #Generate config template
        import string
        config_template = string.Template(open(input.config_template, "r").read())
        config_fill_values = run_configs[wildcards.run]
        config_fill_values['camparee_output'] = "none"# input.camparee_output
        config_fill_values['molecule_files'] = pathlib.Path(params.molecule_dir).resolve()
        config_fill_values['reference_genome'] = pathlib.Path(input.reference_genome).resolve()
        config = config_template.substitute(config_fill_values)
        with open(params.config, "w") as f:
            f.write(config)
        # Run beers
        beers_cmd =  f"run_beers --configfile {params.config} --directory {params.beers_dir} -c 1 --rerun-incomplete"
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
        "/home/thobr/for_tom/BEERS_noselect/data_UMI/samples/WT4_PolyA/deduped.R1.fastq.gz",
        "/home/thobr/for_tom/BEERS_noselect/data_UMI/samples/WT4_PolyA/deduped.R2.fastq.gz",
    output:
        gc_content = "real_data/WT4_PolyA/gc_content.txt"
    params:
        fastq_files = ["/home/thobr/for_tom/BEERS_noselect/data_UMI/samples/WT4_PolyA/deduped.R1.fastq.gz",
            "/home/thobr/for_tom/BEERS_noselect/data_UMI/samples/WT4_PolyA/deduped.R2.fastq.gz"]
    resources:
        mem_mb = 6_000,
    script:
        "scripts/compute_gc_content.py"

rule gather_gc_content:
    input:
        gc_content = expand("data/{run}/gc_content.txt", run=run_configs.keys()),
    output:
        gc_content = "results/gc_content.txt",
        gc_summary = "results/gc_summary.txt"
    run:
        import pandas
        import numpy
        res = []
        for gc_content in input.gc_content:
            d = pandas.read_csv(gc_content, sep="\t")
            res.append(d)
        data = pandas.concat(res)
        data.to_csv(output.gc_content, index=False)

        bins = numpy.linspace(0,1, 101)
        bins[-1] = 1.001 # make it right-inclusive on the last bin
        def get_GC_content(data):
            distribution,_ = numpy.histogram(data.GC_content, bins)
            return pandas.Series(distribution / numpy.sum(distribution), index=bins[:-1])
        GC_summary = data.groupby(["run", "sample"]).apply(get_GC_content)
        GC_summary.to_csv(output.gc_summary, sep="\t")

rule plot_gc_content:
    input:
        gc_summary = "results/gc_summary.txt",
    output:
        directory("results/gc_content/")
    script:
        "scripts/plot_gc_content.py"

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
        bam = "/home/thobr/for_tom/BEERS_noselect/data_UMI/samples/WT4_PolyA/deduped.Aligned.out.sorted.bam",
        gene_ids = "chosen_transcripts.txt",
        gtf = GTF,
        bai = "/home/thobr/for_tom/BEERS_noselect/data_UMI/samples/WT4_PolyA/deduped.Aligned.out.sorted.bam.bai",
    output:
        cov = "real_data/WT4_PolyA/coverage.txt"
    script:
        "scripts/compute_real_coverage.py"

rule compare_real_sim_cov:
    input:
        "data/all_bias/sample1/coverage_summary.txt",
        "real_data/WT4_PolyA/coverage_summary.txt",
    output:
        outdir = directory("results/compare_real_sim_cov/")
    script:
        "scripts/compare_real_sim_cov.py"

rule compute_seq_bias:
    input:
        flag_file = expand("data/{run}/beers/finished_flag",
                        run = run_configs.keys(),),
    output:
        seq_frequencies = "results/seq_bias/seq_frequencies.json",
    script:
        "scripts/compute_seq_bias.py"

rule plot_seq_bias:
    input:
        seq_frequencies = "results/seq_bias/seq_frequencies.json",
    output:
        fwd_frequencies = "results/seq_bias/fwd_seq_frequencies.png",
        rev_frequencies = "results/seq_bias/rev_seq_frequencies.png",
    script:
        "scripts/plot_seq_bias.py"

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

rule make_browser_tracks:
    input:
        bam_files = expand("data/{run}/sample1/BEERS_output.sorted.bam", run = run_configs.keys()),
        bai_files = expand("data/{run}/sample1/BEERS_output.sorted.bam.bai", run = run_configs.keys())
    output:
        "results/browser_tracks/url.txt",
    params:
        track_dir = "results/browser_tracks/",
    run:
        # region to use:
        # chr16:37,868,389-37,888,858
        chrom = "16"
        start = 37868389
        end = 37888858

        track_dir = pathlib.Path(params.track_dir)
        track_dir.mkdir(exist_ok=True)

        temp_dir = track_dir / "temp"
        temp_dir.mkdir(exist_ok=True)


        hub = (f"hub BEERS2_asses\n"
               f"shortLabel BEERS2\n"
               f"longLabel BEERS2 Assess - Apob only - Gregory Grant ITMAT Bioinformatics Lab\n"
               f"genomesFile genomes.txt\n"
               f"email ggrant@pennmedicine.upenn.edu\n"
               #f"descriptionUrl"
               )
        (track_dir / "hub.txt").write_text(hub)

        genomes = (f"genome mm10\n"
                   f"trackDb mm10/trackDb.txt\n"
                   #f"metaTab mm10/tabSeparatedFile.txt"
                   )
        (track_dir / "genomes.txt").write_text(genomes)

        genome_dir = track_dir / "mm10"
        genome_dir.mkdir(exist_ok=True)

        # Make the coverage files
        for run in run_configs.keys():
            track_head = f'track type=wiggle_0 name="{run}" description="{run}" visibility=full autoScale=on color=70,128,83 maxHeightPixels=100:50:20 graphType=bar priority=20'
            cov_head = f"fixedStep chrom=chr{chrom} start={start} step=1"
            cmd = f'( echo "{track_head}"; echo "{cov_head}"; samtools depth -sa -r {chrom}:{start}-{end} data/{run}/sample1/BEERS_output.bam | cut -f3; ) | cat > {genome_dir}/{run}.wig'
            print(cmd)
            shell(cmd)


        #track_db = [(f"track {sample_id}.{direction}\n"
        #             f"bigDataUrl {sample_id}.{direction}.Unique.bw\n"
        #             f"shortLabel {sample_id}.{direction}\n"
        #             f"longLabel {sample_id}.{direction}\n"
        #             f"visibility full\n"
        #             f"autoScale on\n"
        #             f"alwaysZero on\n"
        #             f"type bigWig\n")
        #                 for sample_id in samples
        #                 for direction in ['forward', 'reverse']]
        #(genome_dir / "trackDb.txt").write_text('\n'.join(track_db))

        #print("Uploading to bucket via rclone")
        #COV_URL = f"BEERS2_EXPERIMENTS/SELECTION_FRAG_AND_PCR_RAMP/{wildcards.batch_full}/COV_TRACKS/"
        #rclone_cmd = f"rclone copy {track_dir} {RCLONE_REMOTE}:{BUCKET_NAME}/{COV_URL}"
        #print(rclone_cmd)
        #shell(rclone_cmd)

        #for sample_id in samples:
        #    for direction in ['forward', 'reverse']:
        #        bigWig = pathlib.Path(f"../PROCESSED/{wildcards.batch_full}/NORMALIZED_DATA/EXON_INTRON_JUNCTION/COV/{sample_id}.{direction}.Unique.bw")
        #        rclone_cmd = f"rclone copyto {bigWig} {RCLONE_REMOTE}:{BUCKET_NAME}/{COV_URL}mm10/{bigWig.name}"
        #        print(rclone_cmd)
        #        shell(rclone_cmd)

        ## Write out the URL of the trackhub for easy use in the genome browser
        #(track_dir / "url.txt").write_text(f"https://{BUCKET_NAME}.s3.amazonaws.com/{COV_URL}hub.txt\n")