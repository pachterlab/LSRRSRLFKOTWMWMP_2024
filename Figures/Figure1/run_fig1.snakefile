# Author: Diane Trout
# Port of analysis script by Rebecca Loving

from contextlib import chdir
import gzip
from pathlib import Path
import shutil
import sys
import multiprocessing

genomes = ["casteij", "mm39"]
kmers = [31, 63]

genome_references = "/home/rloving/ref/"

threads = multiprocessing.cpu_count()

orientation=["f", "r", "rc", "c"]
read_types=["barcode", "umi", "cDNA"]

cells_tcc_mtx = ["cells_x_tcc.barcodes.txt",
                 "cells_x_tcc.ec.txt",
                 "cells_x_tcc.mtx"]
abundance_mtx = ["abundance_1/abundance.gene.tsv",
                 "abundance_1/abundance.h5",
                 "abundance_1/abundance.tsv",
                 "genes.txt",
                 "matrix.abundance.gene.mtx",
                 "matrix.abundance.gene.tpm.mtx",
                 "matrix.abundance.mtx",
                 "matrix.abundance.tpm.mtx",
                 "transcripts.txt"]
cells_genes_mtx = ["cells_x_genes.barcodes.txt", "cells_x_genes.genes.txt",
                   "cells_x_genes.genes.names.txt", "cells_x_genes.mtx",]

configfile: "config.yaml"

rule ALL:
    input:
#        "13G_orientation/rc_barcode_reversed.fastq.gz", #md5:bc6fbe835a5d36a7a2d831ceb84d9e47
#        "13G_splitcode/f_barcode.fastq.gz",   #md5:8880f2a8edd9480b9270ab9b4b36681a
#        "13G_orientation/c_barcode.fastq.gz", #md5:078f27eee71b8c9a97b011c1b355990e
#        "13G_orientation/r_barcode.fastq.gz", #md5:(has a barcode with a starting \n)
#        "13G_merged/13G_barcode.fastq.gz",
#        "13G_merged/13G_cDNA.fastq.gz",
#        "13G_merged/13G_umi.fastq.gz",
#
        # Nanopore 13G single cell by primer type
        expand("b01_nanopore_13G_single_cell_k63_both_mm39/counts_unfiltered/{fn}", fn=cells_genes_mtx),
        expand("b01_nanopore_13G_single_cell_k63_polyT_mm39/counts_unfiltered/{fn}", fn=cells_genes_mtx),
        expand("b01_nanopore_13G_single_cell_k63_randO_mm39/counts_unfiltered/{fn}", fn=cells_genes_mtx),

        # nanopore bulk 13G
        # expand("b01_nanopore_13G_bulk_k63_casteij/counts_unfiltered/{fn}", fn=cells_tcc_mtx),
        expand("b01_nanopore_13G_bulk_k63_casteij/quant_unfiltered/{fn}", fn=abundance_mtx),
        # expand("b01_nanopore_13G_bulk_k63_mm39/counts_unfiltered/{fn}", fn=cells_tcc_mtx),
        expand("b01_nanopore_13G_bulk_k63_mm39/quant_unfiltered/{fn}", fn=abundance_mtx),

        # nanopore bulk 13H
        # expand("b01_nanopore_13H_bulk_k63_casteij/counts_unfiltered/{fn}", fn=cells_tcc_mtx),
        expand("b01_nanopore_13H_bulk_k63_casteij/quant_unfiltered/{fn}", fn=abundance_mtx),
        # expand("b01_nanopore_13H_bulk_k63_mm39/counts_unfiltered/{fn}", fn=cells_tcc_mtx),
        expand("b01_nanopore_13H_bulk_k63_mm39/quant_unfiltered/{fn}", fn=abundance_mtx),

        # illumina bulk 13G
        # expand("b01_next1_13G_bulk_k31_casteij/counts_unfiltered/{fn}", fn=cells_tcc_mtx),
        expand("b01_next1_13G_bulk_k31_casteij/quant_unfiltered/{fn}", fn=abundance_mtx),
        # expand("b01_next1_13G_bulk_k31_mm39/counts_unfiltered/{fn}", fn=cells_tcc_mtx),
        expand("b01_next1_13G_bulk_k31_mm39/quant_unfiltered/{fn}", fn=abundance_mtx),

        # illumina bulk 13H
        # expand("b01_next1_13H_bulk_k31_casteij/counts_unfiltered/{fn}", fn=cells_tcc_mtx),
        expand("b01_next1_13H_bulk_k31_casteij/quant_unfiltered/{fn}", fn=abundance_mtx),
        # expand("b01_next1_13H_bulk_k31_mm39/counts_unfiltered/{fn}", fn=cells_tcc_mtx),
        expand("b01_next1_13H_bulk_k31_mm39/quant_unfiltered/{fn}", fn=abundance_mtx),

        # illumina 13G single cell
        expand("b01_next1_13G_single_cell_k31_both_mm39/counts_unfiltered/{fn}", fn=cells_genes_mtx),
        expand("b01_next1_13G_single_cell_k31_polyT_mm39/counts_unfiltered/{fn}", fn=cells_genes_mtx),
        expand("b01_next1_13G_single_cell_k31_randO_mm39/counts_unfiltered/{fn}", fn=cells_genes_mtx),

# stage parse barcode files
rule copy_parse_barcodes:
    input:
        config["parse_r1_R_barcodes"],
        config["parse_r1_T_barcodes"],
        config["parse_r2_3_barcodes"],
    output:
        temp("{target}/r1_R.txt"),
        temp("{target}/r1_T.txt"),
        temp("{target}/r2_3.txt"),
    run:
        target = Path(output[0]).parent
        if not target.exists():
            print("Creating directory {}".format(target))
            target.mkdir()

        for i, o in zip(input, output):
            print("Copying {i} to {o}".format(i=str(i), o=str(o)))
            shutil.copyfile(i, o)

# build references
rule build_ref:
    input:
        fasta = genome_references + "{genome}.fa.gz",
        gtf = genome_references +"{genome}.gtf.gz",
    output:
        cdna = "ref_{genome}_k{kmer}/{genome}_kb_k{kmer}.cdna.fa",
        index = "ref_{genome}_k{kmer}/{genome}_kb_k{kmer}.idx",
        t2g = "ref_{genome}_k{kmer}/{genome}_kb_t2g_k{kmer}.txt",
    params:
        tmp = "ref_{genome}_k{kmer}/tmp"
    conda:
        "envs/kb_python.yaml"
    shell:
        "kb ref -i {output.index} -k {wildcards.kmer} -g {output.t2g} -f1 {output.cdna} {input.fasta} {input.gtf} --opt-off --tmp {params.tmp}"

# Generate seqspec using specified version
rule nanopore_seqspec:
    input:
        "nanopore-seqspec.yaml"
    output:
        "nanopore-splitcode.txt"
    conda:
        "envs/seqspec.yaml"
    shell:
        "seqspec index -m rna -s file -t splitcode  -o {output} {input}"

# Rule for nanopore splitcode
rule run_splitcode_nanopore:
    input:
        fastq = lambda wildcards: config[wildcards.subpool+"_nanopore"],  # Modify as per your input files
        config = Path(rules.nanopore_seqspec.output[0]).absolute(),
    output:
        output_dir = directory("{subpool}_splitcode"),
        fastqs = expand("{{subpool}}_splitcode/{orientation}_{read_type}.fastq.gz",
                        orientation=orientation,
                        read_type=read_types)
    threads: threads
    conda:
        "envs/splitcode.yaml"
    shell:
        "cd {output.output_dir}; splitcode -c {input.config} {input.fastq} --x-only --gzip -t {threads}"

# honestly snakemake might be better for the job of run_rev

rule c_orientation_barcodes:
    input:
        "{subpool}_splitcode/c_barcode.fastq.gz",
    output:
        temp("{subpool}_orientation/c_barcode.fastq.gz"),
    threads: 1
    run:
        with gzip.open(input[0], "rt") as fastqc:
            with gzip.open(output[0], "wt+") as ofastqc:
                for l in fastqc:
                    if l.startswith("@") or l.startswith("+"):
                        ofastqc.write(l)
                    elif l != "\n":
                        ofastqc.write(l[0:8][::-1]+l[8:16][::-1]+l[16:24][::-1]+"\n")
                    else:
                        ofastqc.write(l)

rule r_orientation_barcodes:
    input:
        "{subpool}_splitcode/r_barcode.fastq.gz"
    output:
        temp("{subpool}_orientation/r_barcode.fastq.gz"),
    threads: 1
    run:
        with gzip.open(input[0], "rt") as fastqc:
            with gzip.open(output[0], "wt+") as ofastqc:
                for l in fastqc:
                    if l.startswith("@") or l.startswith("+"):
                        ofastqc.write(l)
                    elif l != "\n":
                        ofastqc.write(l.rstrip()[::-1]+"\n")
                    else:
                        ofastqc.write(l)

rule rc_orientation_barcodes:
    input:
        "{subpool}_splitcode/rc_barcode.fastq.gz"
    output:
        temp("{subpool}_orientation/rc_barcode_reversed.fastq.gz"),
    threads: 1
    run:
        with gzip.open(input[0], "rt") as fastqc:
            with gzip.open(output[0], "wt+") as ofastqc:
                for l in fastqc:
                    if l.startswith("@") or l.startswith("+"):
                        ofastqc.write(l)
                    elif l != "\n":
                        ofastqc.write(l[16:24]+l[8:16]+l[0:8]+"\n")
                    else:
                        ofastqc.write(l)

rule combine_barcode_fastqs:
    input:
        "{subpool}_splitcode/f_barcode.fastq.gz",
        "{subpool}_orientation/c_barcode.fastq.gz",
        "{subpool}_orientation/r_barcode.fastq.gz",
        "{subpool}_orientation/rc_barcode_reversed.fastq.gz",
    output:
        temp("{subpool}_merged/{subpool}_barcode.fastq.gz"),
    shell:
        "cat {input} > {output}"

rule combine_other_fastqs:
    input:
        "{subpool}_splitcode/f_{read_type}.fastq.gz",
        "{subpool}_splitcode/c_{read_type}.fastq.gz",
        "{subpool}_splitcode/r_{read_type}.fastq.gz",
        "{subpool}_splitcode/rc_{read_type}.fastq.gz",
    output:
        "{subpool}_merged/{subpool}_{read_type}.fastq.gz"
    wildcard_constraints:
        read_type="(cDNA)|(umi)"
    shell:
        "cat {input} > {output}"

rule splitcode_correct_nanopore:
    input:
        config = Path("config-correct.txt"),
        fastq = Path("{subpool}_merged/{subpool}_{read_type}.fastq.gz"),
        barcode = Path("{subpool}_merged/{subpool}_barcode.fastq.gz"),
        # these are reference files with parse barcodes
        r1_R = "{subpool}_corrected_{read_type}/r1_R.txt",
        r1_T = "{subpool}_corrected_{read_type}/r1_T.txt",
        r2_3 = "{subpool}_corrected_{read_type}/r2_3.txt",
    output:
        fastq = "{subpool}_corrected_{read_type}/{subpool}_{read_type}.fastq.gz",
        barcode = "{subpool}_corrected_{read_type}/{subpool}_barcode.fastq.gz",
    params:
        output_dir = directory("{subpool}_corrected_{read_type}"),
    wildcard_constraints:
        read_type="(cDNA)|(umi)"
    threads: threads
    conda:
        "envs/splitcode.yaml"
    shell:
        "cd {params.output_dir}; splitcode -c ../{input.config} --nFastqs=2 --select=0 --gzip -o ../{output.fastq} ../{input.fastq} ../{input.barcode} -t {threads} --gzip"

rule splitcode_merge_corrected_nanopore_barcodes:
    input:
        barcode = "{subpool}_corrected_cDNA/{subpool}_barcode.fastq.gz",
    output:
        barcode = "{subpool}_merged_corrected_both/{subpool}_barcode_both.fastq.gz"
    threads: threads
    conda:
        "envs/splitcode.yaml"
    shell:
        "splitcode -c config.mergeRT -o {output.barcode} {input.barcode} -t {threads} --gzip"

rule splitcode_extract_one_nanopore:
    input:
        config = "config_RT_{barcode_type}.txt",
        fastq = "{subpool}_merged/{subpool}_{read_type}.fastq.gz",
        barcode = "{subpool}_merged/{subpool}_barcode.fastq.gz",
        r1_R = "{subpool}_corrected_{read_type}_{barcode_type}/r1_R.txt",
        r1_T = "{subpool}_corrected_{read_type}_{barcode_type}/r1_T.txt",
        r2_3 = "{subpool}_corrected_{read_type}_{barcode_type}/r2_3.txt",
    output:
        fastq = "{subpool}_corrected_{read_type}_{barcode_type}/{subpool}_{read_type}_{barcode_type}.fastq.gz",
        barcode = "{subpool}_corrected_{read_type}_{barcode_type}/{subpool}_barcode_{barcode_type}.fastq.gz",
    params:
        output_dir = directory("{subpool}_corrected_{read_type}_{barcode_type}"),
    wildcard_constraints:
        barcode_type = "(polyT)|(randO)",
        read_type = "(cDNA)|(umi)"
    threads: threads
    conda:
        "envs/splitcode.yaml"
    shell:
        # We're running in a subdirectory because there's some
        # files being written as a side-effect
        "mkdir -p {params.output_dir}; cd {params.output_dir}; splitcode -c ../{input.config} --nFastqs=2 --select=0 --gzip -o ../{output.fastq} ../{input.fastq} ../{input.barcode} -t {threads}"


# Rule for counting with kb cell
rule nanopore_both_barcode_kb_cellxgene_mtx:
    input:
        cDNA = "{subpool}_corrected_cDNA/{subpool}_cDNA.fastq.gz",
        umi = "{subpool}_corrected_umi/{subpool}_umi.fastq.gz",
        barcode = rules.splitcode_merge_corrected_nanopore_barcodes.output.barcode,
        index = rules.build_ref.output.index,
        t2g = rules.build_ref.output.t2g,
    output:
        expand("b01_nanopore_{{subpool}}_single_cell_k{{kmer}}_both_{{genome}}/counts_unfiltered/{fn}", fn=cells_genes_mtx),
        output_dir = directory("b01_nanopore_{subpool}_single_cell_k{kmer}_both_{genome}")
    conda:
        "envs/kb_python.yaml"
    shell:
        "kb count -k {wildcards.kmer} --long --threshold 0.8 -i {input.index} -g {input.t2g} -o {output.output_dir} -x '2,0,24:1,0,10:0,0,0' {input.cDNA} {input.umi} {input.barcode} --opt-off --mm"

# Rule for counting with kb cell
rule nanopore_one_barcode_kb_cellxgene_mtx:
    input:
        cDNA = "{subpool}_corrected_cDNA_{barcode_type}/{subpool}_cDNA_{barcode_type}.fastq.gz",
        umi = "{subpool}_corrected_umi_{barcode_type}/{subpool}_umi_{barcode_type}.fastq.gz",
        #barcode = rules.splitcode_extract_one_nanopore.output.barcode,
        # just pick one, hopefully they're the same
        barcode = "{subpool}_corrected_cDNA_{barcode_type}/{subpool}_barcode_{barcode_type}.fastq.gz",
        index = rules.build_ref.output.index,
        t2g = rules.build_ref.output.t2g,
    output:
        expand("b01_nanopore_{{subpool}}_single_cell_k{{kmer}}_{{barcode_type}}_{{genome}}/counts_unfiltered/{fn}", fn=cells_genes_mtx),
        output_dir = directory("b01_nanopore_{subpool}_single_cell_k{kmer}_{barcode_type}_{genome}")
    wildcard_constraints:
        barcode_type = "(polyT)|(randO)",
        read_type = "(cDNA)|(umi)",
    conda:
        "envs/kb_python.yaml"
    shell:
        "kb count -k {wildcards.kmer} --long --threshold 0.8 -i {input.index} -g {input.t2g} -o {output.output_dir} -x '2,0,24:1,0,10:0,0,0' {input.cDNA} {input.umi} {input.barcode} --opt-off --mm"

rule nanopore_bulk_cell_gene_tcc_mtx:
    input:
        fastq = lambda wildcards: config[f"{wildcards.subpool}_nanopore"],
        index = rules.build_ref.output.index,
        t2g = rules.build_ref.output.t2g,
    output:
        counts = expand("b01_nanopore_{{subpool}}_bulk_k{{kmer}}_{{genome}}/counts_unfiltered/{fn}", fn=cells_tcc_mtx),
        quant = expand("b01_nanopore_{{subpool}}_bulk_k{{kmer}}_{{genome}}/quant_unfiltered/{fn}", fn=abundance_mtx),
        output_dir = directory("b01_nanopore_{subpool}_bulk_k{kmer}_{genome}"),
    conda:
        "envs/kb_python.yaml"
    shell:
        """kb count --overwrite --long --threshold 0.8  -x 'bulk' \
                    -k {wildcards.kmer} \
                    -i {input.index} \
                    -g {input.t2g} \
                    -o {output.output_dir} \
                    {input.fastq} \
                    --opt-off --parity single --tcc --matrix-to-directories"""

# Illumina

# single cell varying the primers
# is this correcting? or just extracting?
rule splitcode_correct_illumina_fastqs:
    input:
        config = Path("config_next1_{barcode_type}.txt").absolute(),
        fastq_r1 = lambda wildcards: config[wildcards.subpool + "_next1"]["R1"],
        fastq_r2 = lambda wildcards: config[wildcards.subpool + "_next1"]["R2"],
        # these are reference files with parse barcodes
        r1_R = "{subpool}_next1_{barcode_type}/r1_R.txt",
        r1_T = "{subpool}_next1_{barcode_type}/r1_T.txt",
        r2_3 = "{subpool}_next1_{barcode_type}/r2_3.txt",
    output:
        fastq_r1 = "{subpool}_next1_{barcode_type}/{subpool}_next1_{barcode_type}_R1.fastq.gz",
        fastq_r2 = "{subpool}_next1_{barcode_type}/{subpool}_next1_{barcode_type}_R2.fastq.gz",
        fastq_umi = "{subpool}_next1_{barcode_type}/{subpool}_next1_{barcode_type}_umi.fastq.gz",
    params:
        output_dir = "{subpool}_next1_{barcode_type}",
    wildcard_constraints:
        barcode_type = "(polyT)|(randO)|(both)",
    threads: threads
    conda:
        "envs/splitcode.yaml"
    shell:
        "mkdir -p {params.output_dir}; cd {params.output_dir}; splitcode -c {input.config} --nFastqs=2 --select=0 --gzip -o ../{output.fastq_r1} {input.fastq_r1} {input.fastq_r2}  -t {threads} --gzip"

rule splitcode_merge_illumina_fastqs:
    input:
        config = "config.mergeRT.next1",
        fastq_r2 = rules.splitcode_correct_illumina_fastqs.output.fastq_r2,
    output:
        fastq_r2 = "{subpool}_next1_merged_{barcode_type}/{subpool}_next1_R2.fastq.gz"
    wildcard_constraints:
        barcode_type = "(both)",
    threads: threads
    conda:
        "envs/splitcode.yaml"
    shell:
        # merge has no side effects
        "splitcode -c {input.config} -o {output.fastq_r2} {input.fastq_r2} -t {threads} --gzip"

rule illumina_single_cell_genes_mtx_both:
    input:
        fastq_r1 = rules.splitcode_correct_illumina_fastqs.output.fastq_r1,
        fastq_r2 = rules.splitcode_merge_illumina_fastqs.output.fastq_r2,
        fastq_umi = rules.splitcode_correct_illumina_fastqs.output.fastq_umi,
        index = rules.build_ref.output.index,
        t2g = rules.build_ref.output.t2g,
    output:
        expand("b01_next1_{{subpool}}_single_cell_k{{kmer}}_{{barcode_type}}_{{genome}}/counts_unfiltered/{fn}", fn=cells_genes_mtx),
    params:
        output_dir = "b01_next1_{subpool}_single_cell_k{kmer}_{barcode_type}_{genome}"
    wildcard_constraints:
        barcode_type = "(both)"
    conda:
        "envs/kb_python.yaml"
    shell:
        """kb count --overwrite -k {wildcards.kmer} \
                    -i {input.index} \
                    -g {input.t2g} \
                    -o {params.output_dir} \
                    -x '1,0,24:2,0,10:0,0,0' \
                    {input.fastq_r1} {input.fastq_r2} {input.fastq_umi} \
                    --opt-off --mm"""

rule illumina_single_cell_genes_mtx_single:
    input:
        fastq_r1 = rules.splitcode_correct_illumina_fastqs.output.fastq_r1,
        fastq_r2 = rules.splitcode_correct_illumina_fastqs.output.fastq_r2,
        fastq_umi = rules.splitcode_correct_illumina_fastqs.output.fastq_umi,
        index = rules.build_ref.output.index,
        t2g = rules.build_ref.output.t2g,
    output:
        expand("b01_next1_{{subpool}}_single_cell_k{{kmer}}_{{barcode_type}}_{{genome}}/counts_unfiltered/{fn}", fn=cells_genes_mtx),
    params:
        output_dir = "b01_next1_{subpool}_single_cell_k{kmer}_{barcode_type}_{genome}"
    wildcard_constraints:
        barcode_type = "(polyT)|(randO)"
    conda:
        "envs/kb_python.yaml"
    shell:
        """kb count --overwrite -k {wildcards.kmer} \
                    -i {input.index} \
                    -g {input.t2g} \
                    -o {params.output_dir} \
                    -x '1,0,24:2,0,10:0,0,0' \
                    {input.fastq_r1} {input.fastq_r2} {input.fastq_umi} \
                    --opt-off --mm"""

rule illumina_bulk_cells_tcc_mtx:
    input:
        fastq_r1 = lambda wildcards: config[wildcards.subpool + "_next1"]["R1"],
        index = rules.build_ref.output.index,
        t2g = rules.build_ref.output.t2g,
    output:
        counts = expand("b01_next1_{{subpool}}_bulk_k{{kmer}}_{{genome}}/counts_unfiltered/{fn}", fn=cells_tcc_mtx),
        quants = expand("b01_next1_{{subpool}}_bulk_k{{kmer}}_{{genome}}/quant_unfiltered/{fn}", fn=abundance_mtx),
    params:
        output_dir = "b01_next1_{subpool}_bulk_k{kmer}_{genome}"
    conda:
        "envs/kb_python.yaml"
    shell:
        """kb count --overwrite -k {wildcards.kmer} -x 'bulk' \
                    -i {input.index} \
                    -g {input.t2g} \
                    -o {params.output_dir} \
                    {input.fastq_r1} \
                    --opt-off --parity single --tcc --matrix-to-directories"""
