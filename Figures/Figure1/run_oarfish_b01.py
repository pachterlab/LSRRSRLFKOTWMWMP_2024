import os

def run(cdna_ref, bs, fastq):
    os.system("./LRGASP_data/minimap2-2.28_x64-linux/minimap2 -t 32 -ax map-ont -N 100 "+cdna_ref+" "+fastq+" > "+bs+".sam")
    os.system("samtools view -bS -o "+bs+".bam "+bs+".sam")
    os.system("oarfish --alignments "+bs+".bam --threads 32 --output "+bs+"_oarfish --model-coverage --filter-group no-filters")

run("ref/mm39.cdna.fa", "b01_nanopore_13H", "igvfb01_13H_lig-ss.fastq.gz")
run("ref/mm39.cdna.fa", "b01_nanopore_13G", "igvfb01_13G-gc_lig-ss_11.fastq.gz")
