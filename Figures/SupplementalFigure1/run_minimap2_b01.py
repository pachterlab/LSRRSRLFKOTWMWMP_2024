import os

def run(genome_ref, bs, fastq):
    os.system("./LRGASP_data/minimap2-2.28_x64-linux/minimap2 -t 32 -ax splice "+genome_ref+" "+fastq+" > "+bs+".sam")
    os.system("samtools view -bS --fast -@ 16 -o "+bs+".bam "+bs+".sam")
  
run("ref/mm39.fa.gz", "b01_nanopore_13H_spliced", "igvfb01_13H_lig-ss.fastq.gz")
