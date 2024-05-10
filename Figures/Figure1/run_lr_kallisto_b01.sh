reads='igvfb01_13H_lig-ss.fastq.gz' #'../8cube/raw_data/igvf_b01/nanopore/demux/13G_demux.fastq.gz'
reads_2='' #'13G_demux.barcodes.fastq'
#'../8cube/raw_data/igvf_b01/nanopore/demux/13G_demux.fastq.gz' #'igvfb01_13G-gc_lig-ss_11.fastq.gz' #'igvfb01_13H_lig-ss.fastq.gz'
path_to_lr_kallisto='kallisto_synced/kallisto/build/src/kallisto'
path_to_bustools='bustools/build/src/bustools'
output='b01_nanopore_13H_casteij' #G_demux' #H.5'
ref='ref/casteij'
tech='bulk' #'1,0,24:1,24,32:0,0,0' #'bulk'
#kb ref --kallisto ${path_to_lr_kallisto} -i ${ref}_k-63.idx -k 63 -f1 ${ref}.cdna.fa -g ${ref}.t2g ${ref}.fa.gz ${ref}.gtf.gz 
#kb ref --kallisto ${path_to_lr_kallisto} --workflow=nac -i ${ref}_k-63.nac.idx -g ${ref}.t2g \
#  -f1 ${ref}.cdna.fa -f2 ${ref}.nascent.fa -c1 ${ref}.cdna.txt -c2 ${ref}.nascent.txt --overwrite -k 63 \
#  ${ref}.fa.gz ${ref}.gtf.gz

${path_to_lr_kallisto} bus --long --threshold 0.8 --error-rate 0.005 -x ${tech} -i ${ref}_k-63.idx -o ${output}_${tech} ${reads} ${reads_2} -t 32
output=${output}_${tech}

${path_to_bustools} sort -t 32 ${output}/output.bus \
 -o ${output}/sorted.bus; \
 ${path_to_bustools} count ${output}/sorted.bus \
 -t ${output}/transcripts.txt \
 -e ${output}/matrix.ec \
 -o ${output}/count --cm -m \
 -g ${ref}.t2g; \
 ${path_to_lr_kallisto} quant-tcc -t 32 --long -P ONT \
 ${output}/count.mtx \
 -i ${ref}_k-63.idx \
 -f ${output}/flens.txt \
 -e ${output}/count.ec.txt \
 -o ${output}; 
