### ONT
#
#splitcode v0.30.0
#kallisto v0.51.0
#bustools v0.43.2
#
#bambu v3.4.1
#IsoQuant v3.3.0
#oarfish v0.5.1
#
splitcode/build/src/splitcode -c config.txt igvfb01_13G-gc_lig-ss_11.fastq.gz -o igvfb01_13G_scmodified/igvfb01_13G_lig-ss.fastq.gz

python3 run_rev.py

splitcode/build/src/splitcode -c config_RT.txt --nFastqs=2 --select=0 --gzip -o 13G_read_rt.fastq.gz 13G_read.fastq 13G_bc.fastq -t 32;

splitcode/build/src/splitcode -c config_RT.txt --nFastqs=2 --select=0 --gzip -o 13G_umi_rt.fastq.gz 13G_umi.fastq 13G_bc.fastq -t 32;

splitcode/build/src/splitcode -c 13H-config_RT.txt --nFastqs=2 --select=0 --gzip -o 13H_read_rt.fastq.gz 13H_read.fastq 13H_bc.fastq -t 32;

splitcode/build/src/splitcode -c 13H-config_RT.txt --nFastqs=2 --select=0 --gzip -o 13H_umi_rt.fastq.gz 13H_umi.fastq 13H_bc.fastq -t 32;

./run_lrk_sc_new.sh 13G_
./run_lrk_sc_new.sh 13H_

splitcode/build/src/splitcode -c config_RT_polyT.txt --nFastqs=2 --select=0 --gzip -o 13G_read_rt_polyT.fastq.gz 13G_read.fastq 13G_bc.fastq -t 32;

splitcode/build/src/splitcode -c config_RT_polyT.txt --nFastqs=2 --select=0 --gzip -o 13G_umi_rt_polyT.fastq.gz 13G_umi.fastq 13G_bc.fastq -t 32;

./run_lrk_sc_new_polyT.sh

splitcode/build/src/splitcode -c config_RT_randO.txt --nFastqs=2 --select=0 --gzip -o 13G_read_rt_randO.fastq.gz 13G_read.fastq 13G_bc.fastq -t 32;

splitcode/build/src/splitcode -c config_RT_randO.txt --nFastqs=2 --select=0 --gzip -o 13G_umi_rt_randO.fastq.gz 13G_umi.fastq 13G_bc.fastq -t 32;

./run_lrk_sc_new_randO.sh
#
##bulk
#
#./run_lrk_bulk.sh file sample ref
./run_lrk_bulk.sh igvfb01_13G-gc_lig-ss_11.fastq.gz 13G mm39
#
./run_lrk_bulk.sh igvfb01_13H_lig-ss.fastq.gz 13H mm39
#
./run_lrk_bulk.sh igvfb01_13G-gc_lig-ss_11.fastq.gz 13G casteij
#
./run_lrk_bulk.sh igvfb01_13G-gc_lig-ss.fastq.gz 13H casteij
#
### Illumina
#
splitcode/build/src/splitcode -c config_next1.txt --nFastqs=2 --select=0 --gzip -o 13G_next1_R1.fastq.gz ../8cube/raw_data/igvf_b01/next1/B01_13G_R1.fastq.gz ../8cube/raw_data/igvf_b01/next1/B01_13G_R2.fastq.gz -t 32; ./run_k_sc_rt.sh
#
./run_k_sc_rt.sh
#
splitcode/build/src/splitcode -c config_next1_randO.txt --nFastqs=2 --select=0 --gzip -o 13G_next1_randO_R1.fastq.gz ../8cube/raw_data/igvf_b01/next1/B01_13G_R1.fastq.gz ../8cube/raw_data/igvf_b01/next1/B01_13G_R2.fastq.gz -t 32; 
#
./run_k_sc_randO.sh
#
splitcode/build/src/splitcode -c config_next1_polyT.txt --nFastqs=2 --select=0 --gzip -o 13G_next1_polyT_R1.fastq.gz ../8cube/raw_data/igvf_b01/next1/B01_13G_R1.fastq.gz ../8cube/raw_data/igvf_b01/next1/B01_13G_R2.fastq.gz -t 32; 
#
./run_k_sc_polyT.sh
#
#bulk
#
#./run_k_bulk.sh sample ref
./run_k_bulk.sh 13G mm39 
#
./run_k_bulk.sh 13G casteij 
#
./run_k_bulk.sh 13H mm39 
#
./run_k_bulk.sh 13H casteij 
