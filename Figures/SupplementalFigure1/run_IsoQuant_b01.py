import os

bam='b01_nanopore_13G_spliced.bam'
outdir='b01_nanopore_13G_IQ'
os.system('samtools sort -@ 16 -o sorted_'+bam+' '+bam)
os.system('samtools index -@ 16 sorted_'+bam)
os.system('isoquant.py --reference ref/mm39.fa.gz --genedb ref/mm39.gtf.gz --complete_genedb --bam sorted_'+ bam +' --data_type nanopore -o '+outdir)
