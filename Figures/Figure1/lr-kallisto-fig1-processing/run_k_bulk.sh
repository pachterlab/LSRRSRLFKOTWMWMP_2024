path_to_lr_kallisto='kallisto_v0.51.0/kallisto/build/src/kallisto'
path_to_bustools='bustools/build/src/bustools'
output='b01_next1_'$1'_'$2 
ref='ref/'$2
tech='bulk'
#kb ref --kallisto ${path_to_lr_kallisto} -i ${ref}_k-63.idx -k 63 -f1 ${ref}.cdna.fa -g ${ref}.t2g ${ref}.fa.gz ${ref}.gtf.gz 
#kb ref --kallisto ${path_to_lr_kallisto} --workflow=nac -i ${ref}_k-63.nac.idx -g ${ref}.t2g \
#  -f1 ${ref}.cdna.fa -f2 ${ref}.nascent.fa -c1 ${ref}.cdna.txt -c2 ${ref}.nascent.txt --overwrite -k 63 \
#  ${ref}.fa.gz ${ref}.gtf.gz

${path_to_lr_kallisto} bus -x ${tech} -i ${ref}_k-31.idx -o ${output}_${tech} -t 32 ../8cube/raw_data/igvf_b01/next1/B01_$1_R1.fastq.gz
output=${output}_${tech}

${path_to_bustools} sort -t 32 ${output}/output.bus \
 -o ${output}/sorted.bus; \
 ${path_to_bustools} whitelist -o ${output}/whitelist.txt \
 ${output}/sorted.bus; \
 ${path_to_bustools} correct -w ${output}/whitelist.txt \
 -o ${output}/corrected.bus ${output}/sorted.bus; \
 ${path_to_bustools} count ${output}/corrected.bus \
 -t ${output}/transcripts.txt \
 -e ${output}/matrix.ec \
 -o ${output}/count -m --cm \
 -g ${ref}.t2g; \
 ${path_to_lr_kallisto} quant-tcc -t 32 \
 ${output}/count.mtx \
 -i ${ref}_k-31.idx \
 -e ${output}/count.ec.txt \
 -o ${output};


${path_to_bustools} count ${output}/corrected.bus \
 -t ${output}/transcripts.txt \
 -e ${output}/matrix.ec \
 -o ${output}/gcount --em --genecounts \
 -g ${ref}.t2g;
