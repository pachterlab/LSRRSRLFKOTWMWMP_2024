reads='../8cube/raw_data/igvf_b01/nanopore/demux/13G_demux.fastq.gz'
reads_2='13G_demux.barcodes.fastq'
#'../8cube/raw_data/igvf_b01/nanopore/demux/13G_demux.fastq.gz' #'igvfb01_13G-gc_lig-ss_11.fastq.gz' #'igvfb01_13H_lig-ss.fastq.gz'
path_to_lr_kallisto='kallisto_v0.51.0/kallisto/build/src/kallisto'
path_to_bustools='bustools/build/src/bustools'
#output='b01_nanopore_'${sample}'_sc_single_cell' #G_demux' #H.5'
ref='ref/mm39'
tech='2,0,24:1,0,10:0,0,0' #'bulk'
#kb ref --kallisto ${path_to_lr_kallisto} -i ${ref}_k-63.idx -k 63 -f1 ${ref}.cdna.fa -g ${ref}.t2g ${ref}.fa.gz ${ref}.gtf.gz 
#kb ref --kallisto ${path_to_lr_kallisto} --workflow=nac -i ${ref}_k-63.nac.idx -g ${ref}.t2g \
#  -f1 ${ref}.cdna.fa -f2 ${ref}.nascent.fa -c1 ${ref}.cdna.txt -c2 ${ref}.nascent.txt --overwrite -k 63 \
#  ${ref}.fa.gz ${ref}.gtf.gz
sample=$1 #'13G_'
output='b01_nanopore_'${sample}'sc_rt_single_cell'

${path_to_lr_kallisto} bus -x ${tech} --long -i ${ref}_k-63.idx -o ${output} -t 32 ${sample}read_rt.fastq.gz ${sample}umi_rt.fastq.gz ${sample}bc_rt.fastq.gz
#output=${output}_${tech}

${path_to_bustools} sort -t 32 ${output}/output.bus \
 -o ${output}/sorted.bus; \
 ${path_to_bustools} whitelist --threshold 200 -o ${output}/whitelist.txt \
 ${output}/sorted.bus; \
 ${path_to_bustools} correct -w ${output}/whitelist.txt \
 -o ${output}/corrected.bus ${output}/sorted.bus; \ 
 ${path_to_bustools} capture -b -c ${output}/whitelist.txt \
 -o ${output}/captured.bus ${output}/corrected.bus; \
 ${path_to_bustools} count ${output}/captured.bus \
 -t ${output}/transcripts.txt \
 -e ${output}/matrix.ec \
 -o ${output}/count -m \
 -g ${ref}.t2g; \
 ${path_to_lr_kallisto} quant-tcc -t 32 \
 ${output}/count.mtx \
 -i ${ref}_k-63.idx \
 -f ${output}/flens.txt \
 -e ${output}/count.ec.txt \
 -o ${output};

${path_to_bustools} count ${output}/corrected.bus \
 -t ${output}/transcripts.txt \
 -e ${output}/matrix.ec \
 -o ${output}/gcount --em --genecounts \
 -g ${ref}.t2g;
