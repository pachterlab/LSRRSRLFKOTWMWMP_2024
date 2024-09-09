import os

samples = ['13H','13G']


for sample in samples:
    files_r = [sample+'frc', sample+'f']
    files_f = [sample+'rr', sample+'r']
    for f in files_r:
        for fastq in ["bc", "read", "umi"]:
            os.system("gunzip "+f+"_"+fastq+".fastq.gz")
        for fastq in ["bc"]: #, "read", "umi"]:
            os.system("python3 rev_umi_bc_reads.py --fastq "+f+"_"+fastq+".fastq --ofastq "+f+"_"+fastq+"_reversed.fastq")

    for fastq in ["bc", "read", "umi"]:
        for f in files_f:
            os.system("gunzip "+f+"_"+fastq+".fastq.gz")
        cmd = "cat "
        for f in files_r: 
            if fastq == "bc":
                cmd += f+'_'+fastq+'_reversed.fastq '
            else:
                cmd += f+'_'+fastq+'.fastq '
        for f in files_f:
            cmd += f+'_'+fastq+'.fastq '
        cmd += '> '+sample+'_'+fastq+'.fastq'
        print(cmd)
        os.system(cmd)
